% function to create CFS stream (byte array) from raw EEG data
% raw data is 4 X Sample matrix, the 4 channels in order are:
% C3, C4, EL and ER respectively
% 
% streamCFS(EEGData, samplingRate, Compress, Hash);
% EEGData is raw 4 channel raw EEG data. Must be sampled at 100Hz or more
% samplingRate is sampling rate in Hz
% Compress is optional argument (default = 1), to control compression
% Hash is optional argument (default = 1), to control transport security
%
% Patents pending (c)-2016 Amiya Patanaik amiyain@gmail.com
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
% THE SOFTWARE.
function stream = streamCFS(EEGData, samplingRate, varargin)

% set defaults for optional inputs
optargs = {1 1};
numvarargs = length(varargin);
if numvarargs > 2,
    error('streamCFS accepts 4 arguments only.')
end
optargs(1:numvarargs) = varargin;
[compress, hash] = optargs{:};

%Order 50 FIR filter
%Basic Settings
SRATE = 100; %Hz
LOWPASS = 45; %Hz
HIGHPASS = 0.3; %Hz
LOWPASSEOG = 12; %Hz
EPOCH = 30*SRATE; %Samples
Fs = samplingRate/2;
bEEG = fir1(50,[HIGHPASS/Fs LOWPASS/Fs]);
bEOG = fir1(50,[HIGHPASS/Fs LOWPASSEOG/Fs]);

eogL = filter(bEOG,1,EEGData(3,:));
eogR = filter(bEOG,1,EEGData(4,:));
eeg = (filter(bEEG,1,EEGData(1,:)) + filter(bEEG,1,EEGData(2,:)))./2;

[p,q] = rat(SRATE/samplingRate);

if(samplingRate ~= 100),
    eogL = resample(eogL,p,q);
    eogR = resample(eogR,p,q);
    eeg =  resample(eeg,p,q);
end

%Extract features from data in each epoch

totalEpochs = floor(size(eeg,2)/30/SRATE);
data = zeros(32,32,3,totalEpochs);


for j=1:totalEpochs,
    
    % FOR EEG-------------------------------------------
    s = spectrogram(eeg(1,(j-1)*EPOCH+1:j*EPOCH),128,38,128,100,'yaxis');
    s(33:end,:) = [];
    data(:,:,1,j) = abs(s);
    
    
    % FOR EOGL-------------------------------------------
    s = spectrogram(eogL(1,(j-1)*EPOCH+1:j*EPOCH),128,38,128,100,'yaxis');
    s(33:end,:) = [];
    data(:,:,2,j) = abs(s);
    
    
    % FOR EOGR-------------------------------------------
    s = spectrogram(eogR(1,(j-1)*EPOCH+1:j*EPOCH),128,38,128,100,'yaxis');
    s(33:end,:) = [];
    data(:,:,3,j) = abs(s);
    
end

%Signature first 3 bytes
signature = uint8('CFS');
%file version next 1 byte
version = uint8(1);
%frequency 1 byte - time 1byte - channel 1byte epochs 2bytes
frequency = uint8(32);
time = uint8(32);
channel = uint8(3);
epochs = typecast(uint16(totalEpochs),'uint8');
%boolean bits
compressionbit = uint8(compress);
hashsetbit = uint8(hash);

%Convert data to single as per specifications
data = cast(data(:),'single')';
actualDigest = [];

if(hashsetbit == 1)
    Opt.Method = 'SHA-1';
    Opt.Format = 'uint8';
    Opt.Input = 'bin';
    actualDigest = DataHash(data, Opt);
end

if(compressionbit == 1)
    data = zlibencode(typecast(data,'uint8'));   
else
    data = typecast(data,'uint8'); 
end

stream = [signature version frequency time channel epochs compressionbit hashsetbit actualDigest data];


end