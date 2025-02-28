#include <iostream>
#include <sndfile.h>
#include <vector>
#include <string>
#include <cstring>
#include <chrono>
#include <cstdlib>
#include <cmath>

#include "logger.hpp"

// bandpass
#define UPPER_FREQUENCY_LIMIT 0.02
#define LOWER_FREQUENCY_LIMIT 0
#define DELTA_F 0.1

// notch
#define NOTCH_N 2
#define NOTCH_f0 2e4

using namespace std;
using namespace std::chrono;

void readWavFile(const std::string &inputFile, std::vector<float> &data, SF_INFO &fileInfo, Logger log)
{
    SNDFILE *inFile = sf_open(inputFile.c_str(), SFM_READ, &fileInfo);
    if (!inFile)
    {
        std::cerr << "Error opening input file: " << sf_strerror(NULL) << std::endl;
        exit(1);
    }

    data.resize(fileInfo.frames * fileInfo.channels);
    sf_count_t numFrames = sf_readf_float(inFile, data.data(), fileInfo.frames);
    if (numFrames != fileInfo.frames)
    {
        std::cerr << "Error reading frames from file." << std::endl;
        sf_close(inFile);
        exit(1);
    }

    sf_close(inFile);
    log.output("Successfully read " + to_string(numFrames) + " frames from " + inputFile);
}

void writeWavFile(const string &outputFile, const vector<float> &data, SF_INFO fileInfo, Logger &log)
{
    sf_count_t originalFrames = fileInfo.frames;
    SNDFILE *outFile = sf_open(outputFile.c_str(), SFM_WRITE, &fileInfo);
    if (!outFile)
    {
        std::cerr << "Error opening output file: " << sf_strerror(NULL) << std::endl;
        exit(1);
    }

    sf_count_t numFrames = sf_writef_float(outFile, data.data(), originalFrames);
    if (numFrames != originalFrames)
    {
        std::cerr << "Error writing frames to file." << std::endl;
        sf_close(outFile);
        exit(1);
    }

    sf_close(outFile);
    log.output("Successfully write " + to_string(numFrames) + " frames to " + outputFile);
}

void applyBandPassFilter(vector<float> &inputData, SF_INFO &fileInfo, Logger log)
{
    const float sampleRate = fileInfo.samplerate;
    for (size_t i = 0; i < inputData.size(); ++i)
    {
        float frequency = sampleRate;
        float filterResponse;
        // cout<< inputData[i]<<" ";
        if (abs(inputData[i]) >= LOWER_FREQUENCY_LIMIT && abs(inputData[i]) <= UPPER_FREQUENCY_LIMIT)
        {
            filterResponse = pow(frequency, 2) / (pow(frequency, 2) + pow(DELTA_F, 2));
        }
        else
        {
            filterResponse = 0.0f; // Outside the passband
        }

        // Apply the filter to the input data
        inputData[i] *= filterResponse;
    }
}

void applyNotchFilter(vector<float> &inputData, SF_INFO &fileInfo, Logger log)
{
    // Calculate the filter response for each frequency
    const float sampleRate = fileInfo.samplerate; // Assume fileInfo is globally available
    for (size_t i = 0; i < inputData.size(); ++i)
    {

        float frequency = sampleRate;
        float filterResponse;
        // cout<< inputData[i]<<" ";
        filterResponse = 1.0f / (1.0f + pow(frequency / NOTCH_f0, 2 * NOTCH_N));
        // Apply the filter to the input data
        inputData[i] *= filterResponse;

        // int o = 1;
        // if(o==1){
        //     cout<<filterResponse<<endl;
        //     o = 0;
        // }
    }
}

void applyFIRFilter(vector<float> &inputData, SF_INFO &fileInfo, Logger log)
{
    // Build the base coefficients vector manually
    // vector<float> coefficients = {
    //     0.005f, -0.001f, 0.001f, 0.002f, -0.0015f, 0.003f, 0.0045f, -0.0005f, 0.0025f, -0.0001f,
    //     0.0018f, -0.0018f, 0.004f, -0.0008f, 0.0033f, 0.0005f, 0.0021f, -0.0012f, 0.0049f, -0.0002f,
    //     0.0017f, -0.0011f, 0.0032f, 0.0008f, 0.0043f, -0.0017f, 0.0027f, -0.0004f, 0.0036f, 0.001f,
    //     0.002f, -0.001f, 0.001f, 0.002f, -0.0015f, 0.003f, 0.0045f, -0.0005f, 0.0025f, -0.0001f,
    //     0.0018f, -0.0018f, 0.004f, -0.0008f, 0.0033f, 0.0005f, 0.0021f, -0.0012f, 0.0049f, -0.0002f,
    //     0.0017f, -0.0011f, 0.0032f, 0.0008f, 0.0043f, -0.0017f, 0.0027f, -0.0004f, 0.0036f, 0.001f,
    //     0.002f, -0.001f, 0.001f, 0.002f, -0.0015f, 0.003f, 0.0045f, -0.0005f, 0.0025f, -0.0001f,
    //     0.0018f, -0.0018f, 0.004f, -0.0008f, 0.0033f, 0.0005f, 0.0021f, -0.0012f, 0.0049f, -0.0002f,
    //     0.0017f, -0.0011f, 0.0032f, 0.0008f, 0.0043f, -0.0017f, 0.0027f, -0.0004f, 0.0036f, 0.001f};
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

    for (size_t i = 0; i < inputData.size(); i++)
    {
        for (size_t k = 0; k < M; k++)
        {
            if (i >= k)
            {
                filteredData[i] += coefficients[k] * inputData[i - k];
            }
        }
    }
    inputData = filteredData;
}

void applyIIRFilter(vector<float> &inputData, SF_INFO &fileInfo, Logger log)
{
    // Build the base coefficients vector manually
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
    float multiplierB = 0.05f;
    vector<float> bCoefficients;
    for (float value : coefficients)
    {
        bCoefficients.push_back(value * multiplierB);
    }

    float multiplierC = 0.4f;
    vector<float> cCoefficients;
    for (float value : coefficients)
    {
        cCoefficients.push_back(value * multiplierC);
    }

    vector<float> filteredData(inputData.size(), 0.0f);
    size_t M = bCoefficients.size();
    size_t N = cCoefficients.size();

    for (size_t i = 0; i < inputData.size(); i++)
    {
        for (size_t k = 0; k < M; k++)
        {
            if (i >= k)
            {
                filteredData[i] += bCoefficients[k] * inputData[i - k];
            }
        }
        // cout << filteredData[i] << " ";

        for (size_t k = 1; k < N; k++)
        {
            if (i >= k)
            {
                filteredData[i] += cCoefficients[k] * filteredData[k - i]*0.1;
            }
        }
    }
        // cout << filteredData[i] << endl;
    inputData = filteredData;
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
    chrono::duration<double> readTime = endRead - startRead;
    log.minfo("Read: " + to_string(readTime.count() * 1000) + " ms");

    // BandPass
    vector<float> bandPassData = audioData;
    auto startBandPass = high_resolution_clock::now();
    applyBandPassFilter(bandPassData, fileInfo, log);
    auto endBandPass = high_resolution_clock::now();
    duration<double> bandPassTime = endBandPass - startBandPass;
    log.ginfo("Bandpass Filter: " + to_string(bandPassTime.count() * 1000) + " ms");
    writeWavFile("outputBandPassSerial.wav", bandPassData, fileInfo, log);

    // Notch
    vector<float> NotchData = audioData;
    auto startNotch = high_resolution_clock::now();
    applyNotchFilter(NotchData, fileInfo, log);
    auto endNotch = high_resolution_clock::now();
    duration<double> NotchTime = endNotch - startNotch;
    log.ginfo("Notch Filter: " + to_string(NotchTime.count() * 1000) + " ms");
    writeWavFile("outputNotchSerial.wav", NotchData, fileInfo, log);

    // FIR
    vector<float> FIRData = audioData;
    auto startFIR = high_resolution_clock::now();
    applyFIRFilter(FIRData, fileInfo, log);
    auto endFIR = high_resolution_clock::now();
    duration<double> FIRTime = endFIR - startFIR;
    log.ginfo("FIR Filter: " + to_string(FIRTime.count() * 1000) + " ms");
    writeWavFile("outputFIRSerial.wav", FIRData, fileInfo, log);

    // IIR
    vector<float> IIRData = audioData;
    auto startIIR = high_resolution_clock::now();
    applyIIRFilter(IIRData, fileInfo, log);
    auto endIIR = high_resolution_clock::now();
    duration<double> IIRTime = endIIR - startIIR;
    log.ginfo("IIR Filter: " + to_string(IIRTime.count() * 1000) + " ms");
    writeWavFile("outputIIRSerial.wav", IIRData, fileInfo, log);

    // final reports
    auto endTotal = high_resolution_clock::now();
    chrono::duration<double> totalTime = endTotal - startTotal;
    log.minfo("Execution: " + to_string(totalTime.count() * 1000) + " ms");

    return 0;
}
