# matcfslib

matcfslib is a MATLAB library to read, write and create compressed feature set (CFS) file/stream from raw PSG data. The CFS format is used by the Z3Score sleep scoring system (https://z3score.com). Instead of using polysomnography data in European Data Format (EDF, https://en.wikipedia.org/wiki/European_Data_Format), the Z3Score system uses CFS files. CFS files are on an average 17X smaller than corresponding EDF files. This reduces data overhead significantly. The format does not allow any user identifiable information ensuring anonymity. The specifications of the file will be made available very soon. The code is released under GPL v3. For alternative license please contact contact@z3score.com 

(c)-2017 Amiya Patanaik 

### Installation

Just add all the folder to Matlab path

### Usage
Z3Score provides a RESTful API to access the sleep scoring services. Read about the API here: https://github.com/amiyapatanaik/z3score-api 
Sample code using the CFS library is provided in the repository. 
### Important Functions

```matlab
    writeCFS(fileName, EEGData, samplingRate, varargin)
```
  - Saves the CFS stream in file_name
  - file_name: is the file name where you want to store the CFS stream. Should have .cfs extension.
  - EEG_data: is a 4 channels X N sample array. The 4 channel in order are C3-A1, C4-A2, EoGleft-A1 and EoGright-A2. Data must be sampled at 100 Hz or more. 
  - sampling_rate: is the signal sampling rate in Hz. All 4 channels must be sampled at the same rate.
  - compressionbit: is True (default) if compression is enabled, False otherwise
  - hashbit: is True (default) if a payload SHA1 signature is included in the CFS stream, False otherwise

```matlab
    stream = streamCFS(EEGData, samplingRate, varargin)
```
  - Returns a CFS binary stream
  - EEG_data: is a 4 channels X N sample array. The 4 channel in order are C3-A1, C4-A2, EoGleft-A1 and EoGright-A2. Data must be sampled at 100 Hz or more. 
  - sampling_rate: is the signal sampling rate in Hz. All 4 channels must be sampled at the same rate.
  - compressionbit: is True (default) if compression is enabled, False otherwise
  - hashbit: is True (default) if a payload SHA1 signature is included in the CFS stream, False otherwise

```matlab
    data = readCFS(filename)
```
  - Returns the data as a 4D array (frequencyXtimeXchannelXepochs) 
  - stream: is the CFS data stream

```matlab
    data = readCFS(filename)
```
  - Returns the data as a 4D array (frequencyXtimeXchannelXepochs) 
  - stream: is the CFS data stream


```matlab
    EDF2CFS(filepath):
```
  - Converts an EDF file to CFS file

 
License
----

GPL V3
