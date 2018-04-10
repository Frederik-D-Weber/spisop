% function to read CFS files to Matlab
% Patents pending (c)-2016 Amiya Patanaik amiyain@gmail.com
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
% THE SOFTWARE.
function data = readCFS(filename)
    data = [];
    fileID = fopen(filename);
    if(fileID == -1),
        disp('File not found!');
        return;
    end
    %Read Header Info
    %Signature first 3 bytes
    signature = fread(fileID,3,'*char','ieee-le');
    if(strcmp(signature,'CFS')~=0),
        disp('File is not a valid CFS.');
        return;
    end
    %file version next 1 byte
    version = fread(fileID,1);
    if(version ~= 1),
        disp('Invalid version number');
        return;
    end
    
    %frequency 1 byte - time 1byte - channel 1byte epochs 2bytes
    freq = fread(fileID,1);
    time = fread(fileID,1);
    channel = fread(fileID,1);
    epochs = fread(fileID,1,'*int16','ieee-le');
    
    %boolean bits
    compressionbit =  fread(fileID,1);
    hashsetbit =  fread(fileID,1);
    
    %Now read the SHA-1 Hash (20 bytes)
    if(hashsetbit == 1)
        digest = fread(fileID,20,'*uint8','ieee-le');
    end
    
    %Now read the raw stream
    stream = uint8(fread(fileID));
    if(compressionbit == 1)
        stream = typecast(zlibdecode(stream),'single');
    else
        stream = typecast(stream,'single');
    end
    
    if(hashsetbit == 1)
        Opt.Method = 'SHA-1';
        Opt.Format = 'uint8';
        Opt.Input = 'bin';
        actualDigest = DataHash(stream, Opt)';
        if(any(actualDigest~=digest)),
            disp('File is corrupt.')
            return;
        end
    end
    
    data = double(reshape(stream,freq,time,channel,epochs));
    fclose(fileID);
end