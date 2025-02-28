#include <iostream>
#include <sndfile.h>
#include <vector>
#include <string>
#include <cstring>
#include <chrono>
#include <cstdlib>
#include <cmath>
#include <thread>
#include <functional>
#include <mutex>

#include "logger.hpp"

// Bandpass Filter 
#define UPPER_FREQUENCY_LIMIT 0.02
#define LOWER_FREQUENCY_LIMIT 0
#define DELTA_F 0.1

// Notch Filter 
#define NOTCH_N 2
#define NOTCH_f0 2e4

using namespace std;
using namespace std::chrono;

// Number of threads allocated for each operation
const size_t READ_THREADS = 12;
const size_t BANDPASS_THREADS = 1;
const size_t NOTCH_THREADS = 1;
const size_t FIR_THREADS = 5;
const size_t IIR_THREADS = 5;

void parallelFor(size_t totalSize, size_t numThreads, function<void(size_t, size_t)> func)
{
    vector<thread> threads;
    size_t chunkSize = totalSize / numThreads;
    for (size_t t = 0; t < numThreads; ++t)
    {
        size_t start = t * chunkSize;
        size_t end = (t == numThreads - 1) ? totalSize : start + chunkSize;
        threads.emplace_back(func, start, end);
    }
    for (auto &th : threads)
    {
        th.join();
    }
}

void readWavFile(const std::string &inputFile, std::vector<float> &data, SF_INFO &fileInfo, Logger &log)
{
    SNDFILE *inFile = sf_open(inputFile.c_str(), SFM_READ, &fileInfo);
    if (!inFile)
    {
        std::cerr << "Error opening input file: " << sf_strerror(NULL) << std::endl;
        exit(1);
    }
    sf_count_t totalFrames = fileInfo.frames;
    sf_close(inFile);

    data.resize(totalFrames * fileInfo.channels, 0.0f);

    // Lambda to read a chunk of data
    auto readChunk = [&](size_t start, size_t end)
    {
        SNDFILE *threadFile = sf_open(inputFile.c_str(), SFM_READ, &fileInfo);
        sf_count_t framesToRead = end - start;
        sf_readf_float(threadFile, data.data() + start * fileInfo.channels, framesToRead);
        sf_close(threadFile);
    };

    parallelFor(totalFrames, READ_THREADS, readChunk);
    log.output("Successfully read " + to_string(totalFrames) + " frames from " + inputFile + " using " + to_string(READ_THREADS) + " threads.");
}

void writeWavFile(const string &outputFile, const vector<float> &data, SF_INFO fileInfo, Logger &log, int Threads_num)
{
    sf_count_t totalFrames = fileInfo.frames;

    SNDFILE *outFile = sf_open(outputFile.c_str(), SFM_WRITE, &fileInfo);
    sf_close(outFile);

    // Lambda to write a chunk of data
    auto writeChunk = [&](size_t start, size_t end)
    {
        SNDFILE *threadFile = sf_open(outputFile.c_str(), SFM_WRITE, &fileInfo);
        // Write frames
        sf_count_t framesToWrite = end - start;
        sf_writef_float(threadFile, data.data() + start * fileInfo.channels, framesToWrite);
        sf_close(threadFile);
    };
    parallelFor(totalFrames, Threads_num, writeChunk);
    log.output("Successfully wrote " + to_string(totalFrames) + " frames to " + outputFile + " using " + to_string(Threads_num) + " threads.");
}

void applyBandPassFilter(vector<float> &inputData, SF_INFO &fileInfo, Logger &log)
{
    const float sampleRate = fileInfo.samplerate;
    size_t totalSize = inputData.size();

    parallelFor(totalSize, BANDPASS_THREADS, [&](size_t start, size_t end)
                {
        for (size_t i = start; i < end; ++i) {
            float frequency = sampleRate;
            float filterResponse;
            if (abs(inputData[i]) >= LOWER_FREQUENCY_LIMIT && abs(inputData[i]) <= UPPER_FREQUENCY_LIMIT) {
                filterResponse = pow(frequency, 2) / (pow(frequency, 2) + pow(DELTA_F, 2));
            } else {
                filterResponse = 0.0f; // Outside the passband
            }

            // Apply the filter to the input data
            inputData[i] *= filterResponse;
        } });

    log.output("Bandpass filter applied using " + to_string(BANDPASS_THREADS) + " threads.");
}

void applyNotchFilter(vector<float> &inputData, SF_INFO &fileInfo, Logger &log)
{
    const float sampleRate = fileInfo.samplerate;
    size_t totalSize = inputData.size();

    parallelFor(totalSize, NOTCH_THREADS, [&](size_t start, size_t end)
                {
        for (size_t i = start; i < end; ++i) {
            float frequency = sampleRate;
            float filterResponse;

            filterResponse = 1.0f / (1.0f + pow(frequency / NOTCH_f0, 2 * NOTCH_N));
            // Apply the filter to the input data
            inputData[i] *= filterResponse;
        } });

    log.output("Notch filter applied using " + to_string(NOTCH_THREADS) + " threads.");
}

void applyFIRFilter(vector<float> &inputData, SF_INFO &fileInfo, Logger &log)
{
    // FIR filter coefficients
    vector<float> coefficients = {
        0.12f, 4.5f, 3.3f, 1.5f, 0.9f, 2.8f, 4.1f, 3.7f, 2.2f, 1.4f,
        0.6f, 5.0f, 4.3f, 0.8f, 2.9f, 3.5f, 1.2f, 0.4f, 2.7f, 1.6f,
        4.8f, 2.5f, 3.9f, 0.7f, 0.0f, 2.3f, 1.8f, 0.2f, 3.4f, 4.7f,
        1.9f, 2.6f, 0.3f, 3.8f, 1.1f, 4.0f, 2.1f, 3.2f, 1.3f, 0.5f,
        4.4f, 2.4f, 0.0f, 0.9f, 4.6f, 1.7f, 0.8f, 3.6f, 2.0f, 4.9f,
        1.4f, 0.7f, 2.9f, 1.3f, 4.2f, 3.3f, 2.6f, 0.6f, 3.8f, 1.0f,
        4.7f, 2.5f, 0.4f, 3.9f, 1.2f, 4.0f, 2.7f, 1.8f, 0.3f, 4.4f,
        2.2f, 3.4f, 1.5f, 4.8f, 0.1f, 3.1f, 2.0f, 4.9f, 1.6f, 3.7f,
        0.2f, 4.6f, 2.3f, 3.5f, 1.7f, 4.3f, 0.5f, 3.0f, 2.1f, 4.1f,
        1.4f, 0.9f, 2.8f, 3.6f, 4.5f, 1.3f, 3.9f, 2.4f, 0.8f, 4.7f};
    vector<float> filteredData(inputData.size(), 0.0f);
    size_t M = coefficients.size();

    parallelFor(inputData.size(), FIR_THREADS, [&](size_t start, size_t end)
                {
        for (size_t i = start; i < end; i++) {
            for (size_t k = 0; k < M; k++) {
                if (i >= k) {
                    filteredData[i] += coefficients[k] * inputData[i - k];
                }
            }
        } });

    inputData = filteredData;
    log.output("FIR filter applied using " + to_string(FIR_THREADS) + " threads.");
}

void applyIIRFilter(vector<float> &inputData, SF_INFO &fileInfo, Logger &log)
{
    // IIR filter coefficients
    vector<float> coefficients = {
        0.12f, 4.5f, 3.3f, 1.5f, 0.9f, 2.8f, 4.1f, 3.7f, 2.2f, 1.4f,
        0.6f, 5.0f, 4.3f, 0.8f, 2.9f, 3.5f, 1.2f, 0.4f, 2.7f, 1.6f,
        4.8f, 2.5f, 3.9f, 0.7f, 0.0f, 2.3f, 1.8f, 0.2f, 3.4f, 4.7f,
        1.9f, 2.6f, 0.3f, 3.8f, 1.1f, 4.0f, 2.1f, 3.2f, 1.3f, 0.5f,
        4.4f, 2.4f, 0.0f, 0.9f, 4.6f, 1.7f, 0.8f, 3.6f, 2.0f, 4.9f,
        1.4f, 0.7f, 2.9f, 1.3f, 4.2f, 3.3f, 2.6f, 0.6f, 3.8f, 1.0f,
        4.7f, 2.5f, 0.4f, 3.9f, 1.2f, 4.0f, 2.7f, 1.8f, 0.3f, 4.4f,
        2.2f, 3.4f, 1.5f, 4.8f, 0.1f, 3.1f, 2.0f, 4.9f, 1.6f, 3.7f,
        0.2f, 4.6f, 2.3f, 3.5f, 1.7f, 4.3f, 0.5f, 3.0f, 2.1f, 4.1f,
        1.4f, 0.9f, 2.8f, 3.6f, 4.5f, 1.3f, 3.9f, 2.4f, 0.8f, 4.7f};
    // float multiplierB = 0.05f;
    // vector<float> bCoefficients;
    // for (float value : coefficients)
    // {
    //     bCoefficients.push_back(value * multiplierB);
    // }

    float multiplierC = 0.4f;
    vector<float> cCoefficients;
    for (float value : coefficients)
    {
        cCoefficients.push_back(value * multiplierC);
    }

    vector<float> filteredData(inputData.size(), 0.0f);
    // size_t M = bCoefficients.size();
    size_t N = cCoefficients.size();


    applyFIRFilter(inputData, fileInfo, log);

    // IIR filter is inherently sequential due to feedback
    for (size_t i = 0; i < inputData.size(); i++)
    {
        for (size_t k = 1; k < N; k++)
        {
            if (i >= k)
            {
                filteredData[i] += cCoefficients[k] * filteredData[i - k] * 0.1f;
            }
        }
    }

    inputData = filteredData;

    log.output("FIR filter applied using " + to_string(IIR_THREADS) + " threads.");
}

int main(int argc, char **argv)
{
    Logger log("MAIN");
    string inputFile = argv[1];

    SF_INFO fileInfo;
    memset(&fileInfo, 0, sizeof(fileInfo));
    auto startTotal = high_resolution_clock::now();

    auto startRead = high_resolution_clock::now();
    vector<float> audioData;
    readWavFile(inputFile, audioData, fileInfo, log);
    auto endRead = high_resolution_clock::now();
    duration<double> readTime = endRead - startRead;
    log.minfo("Read: " + to_string(readTime.count() * 1000) + " ms");

    // BandPass Filter
    vector<float> bandPassData = audioData;
    auto startBandPass = high_resolution_clock::now();
    applyBandPassFilter(bandPassData, fileInfo, log);
    auto endBandPass = high_resolution_clock::now();
    duration<double> bandPassTime = endBandPass - startBandPass;
    log.ginfo("Bandpass Filter: " + to_string(bandPassTime.count() * 1000) + " ms");

    // Writing BandPass Output
    auto startWriteBP = high_resolution_clock::now();
    writeWavFile("outputBandPassParallel.wav", bandPassData, fileInfo, log, BANDPASS_THREADS);
    auto endWriteBP = high_resolution_clock::now();
    duration<double> writeBPTime = endWriteBP - startWriteBP;
    log.minfo("Write BandPass: " + to_string(writeBPTime.count() * 1000) + " ms");

    // Notch Filter
    vector<float> notchData = audioData;
    auto startNotch = high_resolution_clock::now();
    applyNotchFilter(notchData, fileInfo, log);
    auto endNotch = high_resolution_clock::now();
    duration<double> notchTime = endNotch - startNotch;
    log.ginfo("Notch Filter: " + to_string(notchTime.count() * 1000) + " ms");

    // Writing Notch Output
    auto startWriteNotch = high_resolution_clock::now();
    writeWavFile("outputNotchParallel.wav", notchData, fileInfo, log, NOTCH_THREADS);
    auto endWriteNotch = high_resolution_clock::now();
    duration<double> writeNotchTime = endWriteNotch - startWriteNotch;
    log.minfo("Write Notch: " + to_string(writeNotchTime.count() * 1000) + " ms");

    // FIR Filter
    vector<float> FIRData = audioData;
    auto startFIR = high_resolution_clock::now();
    applyFIRFilter(FIRData, fileInfo, log);
    auto endFIR = high_resolution_clock::now();
    duration<double> FIRTime = endFIR - startFIR;
    log.ginfo("FIR Filter: " + to_string(FIRTime.count() * 1000) + " ms");

    // Writing FIR Output
    auto startWriteFIR = high_resolution_clock::now();
    writeWavFile("outputFIRParallel.wav", FIRData, fileInfo, log, FIR_THREADS);
    auto endWriteFIR = high_resolution_clock::now();
    duration<double> writeFIRTime = endWriteFIR - startWriteFIR;
    log.minfo("Write FIR: " + to_string(writeFIRTime.count() * 1000) + " ms");

    // IIR Filter
    vector<float> IIRData = audioData;
    auto startIIR = high_resolution_clock::now();
    applyIIRFilter(IIRData, fileInfo, log);
    auto endIIR = high_resolution_clock::now();
    duration<double> IIRTime = endIIR - startIIR;
    log.ginfo("IIR Filter: " + to_string(IIRTime.count() * 1000) + " ms");

    // Writing IIR Output
    auto startWriteIIR = high_resolution_clock::now();
    writeWavFile("outputIIRParallel.wav", IIRData, fileInfo, log, IIR_THREADS);
    auto endWriteIIR = high_resolution_clock::now();
    duration<double> writeIIRTime = endWriteIIR - startWriteIIR;
    log.minfo("Write IIR: " + to_string(writeIIRTime.count() * 1000) + " ms");

    // Final Reports
    auto endTotal = high_resolution_clock::now();
    duration<double> totalTime = endTotal - startTotal;
    log.minfo("Total Execution: " + to_string(totalTime.count() * 1000) + " ms");

    return 0;
}
