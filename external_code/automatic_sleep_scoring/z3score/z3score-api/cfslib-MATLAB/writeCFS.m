% function to write CFS files from raw EEG data
% raw data is 4 X Sample matrix, the 4 channels in order are:
% C3, C4, EL and ER respectively
% Patents pending (c)-2016 Amiya Patanaik amiyain@gmail.com
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
% THE SOFTWARE.
function writeCFS(fileName, EEGData, samplingRate, varargin)

% set defaults for optional inputs
optargs = {1 1};
numvarargs = length(varargin);
if numvarargs > 2,
    error('streamCFS accepts 4 arguments only.')
end
optargs(1:numvarargs) = varargin;
[compress, hash] = optargs{:};

stream = streamCFS(EEGData, samplingRate, compress, hash);
fileID = fopen(fileName,'w');

if(fileID == -1),
    disp('File not found!');
    return;
end

fwrite(fileID,stream,'*uint8','ieee-le');
fclose(fileID);
end