% Function to connect with the Z3Score XML-RPC server
% Patents pending (c)-2016 Amiya Patanaik amiyain@gmail.com
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
% THE SOFTWARE.

function client = connectZ3Score(server)

if(nargin < 1),
    server = 'http://localhost:18778';
end

client = javaObject('org.apache.xmlrpc.client.XmlRpcClient');
config = javaObject('org.apache.xmlrpc.client.XmlRpcClientConfigImpl');
url = javaObject('java.net.URL',server);
config.setServerURL(url);
client.setConfig(config);

end

