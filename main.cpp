#include "tcp/hdr/TcpServer.h"
#include "tcp/hdr/TcpClient.h"

#include <iostream>
#include <mutex>
int cntr = 0;
//Parse ip to std::string
std::string getHostStr(const TcpServer::Client& client) {
    uint32_t ip = client.getHost ();
    return std::string() + std::to_string(int(reinterpret_cast<char*>(&ip)[0])) + '.' +
            std::to_string(int(reinterpret_cast<char*>(&ip)[1])) + '.' +
            std::to_string(int(reinterpret_cast<char*>(&ip)[2])) + '.' +
            std::to_string(int(reinterpret_cast<char*>(&ip)[3])) + ':' +
            std::to_string( client.getPort ());
}

TcpServer server( 4002,

[](DataBuffer data, TcpServer::Client& client){ // Data handler
  std::cout << cntr<<" Client "<<getHostStr(client)<<" send data [ " << data.size << " bytes ]: "/* << (char*)data.data_ptr */;
  float tmp = *(float*)(data.data_ptr);
  std::cout << tmp<<",  ";
  for(int i=0;i<10;i++)
  {
    tmp = *(float*)(data.data_ptr+i*sizeof(float));
   std::cout << tmp<<",  ";
  }
  std::cout << tmp<<'\n';

      /*std::cout.unsetf(std::ios::hex);
    std::cout.setf(std::ios::dec);*/
  cntr++;
  client.sendData("Hello, client!", sizeof("Hello, client!"));
},

[](TcpServer::Client& client) { // Connect handler
  std::cout << "Client " << getHostStr(client) << " connected\n";
},


[](TcpServer::Client& client) { // Disconnect handler
  std::cout << "Client " << getHostStr(client) << " disconnected\n";
},

{1, 1, 1} // Keep alive{idle:1s, interval: 1s, pk_count: 1}
);

void testServer() {
  //Start server
  if(server.start() == TcpServer::status::up) {
      std::cout<<"Server listen on port:"<<server.getPort()<<std::endl;
      server.joinLoop();
  } else {
      std::cout<<"Server start error! Error code:"<< int(server.getStatus()) <<std::endl;
  }
}




int main() {
  using namespace std::chrono_literals;
  try {
  testServer();

  std::this_thread::sleep_for(30s);
  } catch(std::exception& except) {
    std::cerr << except.what();
  }
}
