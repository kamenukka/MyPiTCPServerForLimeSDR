#include "tcp/hdr/TcpServer.h"
#include "tcp/hdr/TcpClient.h"

#include <iostream>
#include <mutex>

#include <fstream>

#include <chrono>
#include <math.h>
#include "lime/LimeSuite.h"
#include "myFunc/funcForStudents.h"
#include "myFunc/complex.h"

#define TRANSFER_BLOCK 1024

int cntr = 0;
float tmp,Ku,f1,f2,f3,f4;
float lenArray=0;
int cntrFloatSymbols;
bool flag_start=false;
float rec_buffer[255999996];
int16_t out_buffer[255999996];



float Mode;
float Freq;
float Fd;
float Gain;
float OverSampling;
float NumChanTrans;
float NumChanRec;
float SecondsToPlay;
float SizeOfRecArr;
float tx_streamfifoSize;
float tx_streamthroughputVsLatency;
float rx_streamfifoSize;
float rx_streamthroughputVsLatency;
float fmtData;




std::string strAnswer;
bool flagEndOfProc = true;
int statusOfSDR = 0; 
const complex complex::j(0., 1.);
//Device structure, should be initialize to NULL
lms_device_t* device = NULL;

int error()
{
    if (device != NULL)
  	LMS_Close(device);
    exit(-1);
}

int playDeviceTx()
{
    const double frequency = Freq*1e6;  //center frequency to 500 MHz
    const double sample_rate = Fd*1e3;    //sample rate to 5 MHz
    const double tone_freq = 1.25e4; //tone frequency
    const double f_ratio = tone_freq/sample_rate;
    //Find devices
    int n;
    lms_info_str_t list[8]; //should be large enough to hold all detected devices
    if ((n = LMS_GetDeviceList(list)) < 0) //NULL can be passed to only get number of devices
        error();

    cout << "Devices found: " << n << endl; //print number of devices
    if (n < 1)
        return -1;

    //open the first device
    if (LMS_Open(&device, list[0], NULL))
        error();

    //Initialize device with default configuration
    //Do not use if you want to keep existing configuration
    //Use LMS_LoadConfig(device, "/path/to/file.ini") to load config from INI
    if (LMS_Init(device)!=0)
        error();

    //Enable TX channel,Channels are numbered starting at 0
    if (LMS_EnableChannel(device, LMS_CH_TX, 0, true)!=0)
        error();

    //Set sample rate
    if (LMS_SetSampleRate(device, sample_rate, 0)!=0)
        error();
    cout << "Sample rate: " << sample_rate/1e6 << " MHz" << endl;

    //Set center frequency
    if (LMS_SetLOFrequency(device,LMS_CH_TX, 0, frequency)!=0)
        error();
    cout << "Center frequency: " << frequency/1e6 << " MHz" << endl;

    //select TX1_1 antenna
    if (LMS_SetAntenna(device, LMS_CH_TX, 0, LMS_PATH_TX1)!=0)
        error();

    //set TX gain
    if (LMS_SetNormalizedGain(device, LMS_CH_TX, 0, Gain) != 0)
        error();

    //calibrate Tx, continue on failure
    LMS_Calibrate(device, LMS_CH_TX, 0, sample_rate, 0);
    
    //Streaming Setup
    
    lms_stream_t tx_stream;                 //stream structure
    tx_stream.channel = 0;                  //channel number
    tx_stream.fifoSize = int(tx_streamfifoSize);          //fifo size in samples
    tx_stream.throughputVsLatency = 0.5;    //0 min latency, 1 max throughput
    tx_stream.dataFmt = lms_stream_t::LMS_FMT_F32; //floating point samples
    tx_stream.isTx = true;                  //TX channel
    LMS_SetupStream(device, &tx_stream);
    int lenFile=cntrFloatSymbols;
    string line;
    
    //Initialize data buffers
    //const int buffer_size = ((cntrFloatSymbols/2)%1024)?((int((cntrFloatSymbols/2)/1024)+1)*1024):int((cntrFloatSymbols/2));
    const int buffer_size = lenFile/2;

    float tx_buffer[2*buffer_size];     //buffer to hold complex values (2*samples))
    for (int i = 0; i <buffer_size*2; i++) {      //generate TX tone
        tx_buffer[i] = 0;
    }  
    for (int i = 0; i <lenFile; i++) {      //generate TX tone
        tx_buffer[i] = rec_buffer[i];
        //cout << tx_buffer[i]<<std::endl;
    }   
    const int send_cnt = int(buffer_size*f_ratio) / f_ratio; 
    cout << "sample count per send call: " << send_cnt << std::endl;

    LMS_StartStream(&tx_stream);         //Start streaming
    //Streaming
    auto t1 = chrono::high_resolution_clock::now();
    auto t2 = t1;
    while (chrono::high_resolution_clock::now() - t1 < chrono::seconds(int(SecondsToPlay))) //run for 10 seconds
    {
        //Transmit samples
        int ret = LMS_SendStream(&tx_stream, tx_buffer, send_cnt, nullptr, 1000);
        if (ret != send_cnt)
            cout << "error: samples sent: " << ret << "/" << send_cnt << endl;
        //Print data rate (once per second)
        if (chrono::high_resolution_clock::now() - t2 > chrono::seconds(1))
        {
            t2 = chrono::high_resolution_clock::now();
            lms_stream_status_t status;
            LMS_GetStreamStatus(&tx_stream, &status);  //Get stream status
            cout << "TX data rate: " << status.linkRate / 1e6 << " MB/s\n"; //link data rate
        }
    }
    //Stop streaming
    LMS_StopStream(&tx_stream);
    LMS_DestroyStream(device, &tx_stream);

    //Disable TX channel
    if (LMS_EnableChannel(device, LMS_CH_TX, 0, false)!=0)
        error();

    //Close device
    if (LMS_Close(device)==0)
        cout << "Closed" << endl;
    return 0;
  }



int recordDeviceRx()
{
    //Find devices
    const double frequency = Freq*1e6;  //center frequency to 500 MHz
    const double sample_rate = Fd*1e3;    //sample rate to 5 MHz
    int n;
    lms_info_str_t list[8]; //should be large enough to hold all detected devices
    if ((n = LMS_GetDeviceList(list)) < 0) //NULL can be passed to only get number of devices
        error();

    cout << "Devices found: " << n << endl; //print number of devices
    if (n < 1)
        return -1;

    //open the first device
    if (LMS_Open(&device, list[0], NULL))
        error();

    //Initialize device with default configuration
    //Do not use if you want to keep existing configuration
    //Use LMS_LoadConfig(device, "/path/to/file.ini") to load config from INI
    if (LMS_Init(device) != 0)
        error();

    //Enable RX channel
    //Channels are numbered starting at 0
    if (LMS_EnableChannel(device, LMS_CH_RX, 0, true) != 0)
        error();

    //Set center frequency to 800 MHz
    if (LMS_SetLOFrequency(device, LMS_CH_RX, 0, frequency) != 0)
        error();

    //Set sample rate to 8 MHz, ask to use 2x oversampling in RF
    //This set sampling rate for all channels
    if (LMS_SetSampleRate(device, sample_rate, OverSampling) != 0)
        error();

    //Enable test signal generation
    //To receive data from RF, remove this line or change signal to LMS_TESTSIG_NONE
    //if (LMS_SetTestSignal(device, LMS_CH_RX, 0, LMS_TESTSIG_NCODIV8, 0, 0) != 0)
        //error();

    //Streaming Setup

    //Initialize stream
    lms_stream_t streamId; //stream structure
    streamId.channel = 0; //channel number
    streamId.fifoSize = 1024 * 1024; //fifo size in samples
    streamId.throughputVsLatency = 1.0; //optimize for max throughput
    streamId.isTx = false; //RX channel
    streamId.dataFmt = lms_stream_t::LMS_FMT_I12; //12-bit integers
    if (LMS_SetupStream(device, &streamId) != 0)
        error();

    //Initialize data buffers
    const int sampleCnt = int(SizeOfRecArr); //complex samples per buffer
    int16_t buffer[sampleCnt * 2]; //buffer to hold complex values (2*samples))

    //Start streaming
    LMS_StartStream(&streamId);

    //Streaming
#ifdef USE_GNU_PLOT
    GNUPlotPipe gp;
#endif
    //int16_t * in;
    //int16_t * out;
    //auto t1 = chrono::high_resolution_clock::now();
    int cnt =0;
    while (cnt++<1) //run for 5 seconds TODO: надо логично сделать
    {
        //Receive samples
	int samplesRead = LMS_RecvStream(&streamId, buffer, sampleCnt, NULL, 1000);
	if (samplesRead == sampleCnt)
	{
        for (int i = 0; i <samplesRead*2; i++)       //generate TX tone
		{
		  out_buffer[i] =  buffer[i];
            //cout << ": " << out_buffer[i] << endl; //print number of devices
		}
	}
    }
    
    //Stop streaming
    LMS_StopStream(&streamId); //stream is stopped but can be started again with LMS_StartStream()
    LMS_DestroyStream(device, &streamId); //stream is deallocated and can no longer be used

    //Close device
    LMS_Close(device);

    return 0;
}




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
  //std::cout << cntr<<" Client "<<getHostStr(client)<<" send data [ " << data.size << " bytes ]: "/* << (char*)data.data_ptr */;

   if ((*(float*)(data.data_ptr+0*sizeof(float))==1.0)&&
      (*(float*)(data.data_ptr+1*sizeof(float))==1.0)&&
      (*(float*)(data.data_ptr+2*sizeof(float))==0.0)&&
      (*(float*)(data.data_ptr+3*sizeof(float))==0.0)&&
      (*(float*)(data.data_ptr+4*sizeof(float))==1.0)&&
      (*(float*)(data.data_ptr+5*sizeof(float))==1.0))
  {
    flag_start = true;
    flagEndOfProc = false;
    cntr = 0;
    cntrFloatSymbols = 0;
    //place to settings
    Mode = *(float*)(data.data_ptr+ 6*sizeof(float));
    Freq = *(float*)(data.data_ptr+ 8*sizeof(float));
    Fd = *(float*)(data.data_ptr+ 9*sizeof(float));
    Gain = *(float*)(data.data_ptr+ 10*sizeof(float));
    OverSampling = *(float*)(data.data_ptr+ 11*sizeof(float));
    NumChanTrans = *(float*)(data.data_ptr+ 12*sizeof(float));
    NumChanRec = *(float*)(data.data_ptr+ 13*sizeof(float));
    SecondsToPlay = *(float*)(data.data_ptr+ 14*sizeof(float));
    SizeOfRecArr = *(float*)(data.data_ptr+ 15*sizeof(float));
    tx_streamfifoSize = *(float*)(data.data_ptr+ 16*sizeof(float));
    tx_streamthroughputVsLatency = *(float*)(data.data_ptr+ 17*sizeof(float));
    rx_streamfifoSize = *(float*)(data.data_ptr+ 18*sizeof(float));
    rx_streamthroughputVsLatency = *(float*)(data.data_ptr+ 19*sizeof(float));
    fmtData = *(float*)(data.data_ptr+ 20*sizeof(float));
    
    tmp =*(float*)(data.data_ptr+7*sizeof(float));
    lenArray = (tmp);
    std::cout <<"start packet, all length: "<< lenArray<<std::endl;;
    if (Mode)
    {
    for (int k=21;k<(1024/4);k++)
    {
      if (cntrFloatSymbols<lenArray)
        rec_buffer[cntrFloatSymbols++] = *(float*)(data.data_ptr+k*sizeof(float));
        //std::cout<<"!"<<(cntrFloatSymbols++)<<":   "<<*(float*)(data.data_ptr+k*sizeof(float))<<'\n';
      }
     }
     else
     {
        uint8_t      bytes[sizeof(int16_t)];
        int statusOfSDRrec = recordDeviceRx();
        strAnswer = "End recording with status";
        strAnswer.append(std::to_string(statusOfSDRrec));
        const char* tmp2 = strAnswer.c_str();
        bool statSend = client.sendData(tmp2, sizeof(strAnswer)*strAnswer.length());


        int sizeOfTransferBuffer = (SizeOfRecArr)*2*sizeof(int16_t);
        //cout<<sizeof(int16_t);
        char * arrayToTransfer = new char [sizeOfTransferBuffer];
        for (int j=0;j<int(SizeOfRecArr);j++)
        {
          *(int16_t*)(bytes) = out_buffer[j];  // convert float to bytes
          for (int k=0;k<sizeof(int16_t);k++)
              arrayToTransfer[j*sizeof(int16_t)+k]=bytes[k];
         }
        
        std::cout <<"sizeOfTransferBuffer" << sizeOfTransferBuffer << '\n';
        cntr = ((sizeOfTransferBuffer%TRANSFER_BLOCK)!=0)?(int(sizeOfTransferBuffer/TRANSFER_BLOCK)+1):(int(sizeOfTransferBuffer/TRANSFER_BLOCK));
        if (statusOfSDRrec<0)
               cntr = 0;  
        strAnswer = "Server Up with send Data ";
        strAnswer.append(std::to_string(cntr));
                strAnswer.append(" bloks");

        const char* tmp2_ = strAnswer.c_str();
        statSend = client.sendData(tmp2_, sizeof(strAnswer)*strAnswer.length());
        if (statusOfSDRrec<0)
        {
                   strAnswer = "Please check you SDR. Server Ready To Next Command";
    const char* tmp2_ = strAnswer.c_str();
    statSend = client.sendData(tmp2_, sizeof(strAnswer)*strAnswer.length());
           }
           else
           {
        for (int ij=0;ij<cntr;ij++)
        {
         std::cout <<"send " << ij << '\n';
          client.sendData(arrayToTransfer +TRANSFER_BLOCK*ij, TRANSFER_BLOCK);
        }
    }

     }
  }
  else
  {
     /*std::cout <<*(float*)(data.data_ptr+1*sizeof(float))<<", "
               <<*(float*)(data.data_ptr+2*sizeof(float))<<", "
               <<*(float*)(data.data_ptr+3*sizeof(float))<<", "
               <<*(float*)(data.data_ptr+4*sizeof(float))<<", "
               <<*(float*)(data.data_ptr+5*sizeof(float))<<", "
               <<*(float*)(data.data_ptr+6*sizeof(float))<<", "
               <<*(float*)(data.data_ptr+7*sizeof(float))<<", "
               <<*(float*)(data.data_ptr+8*sizeof(float))<<", "
               <<'\n';*/
    for (int k=0;k<(1024/4);k++)
    {
        if (cntrFloatSymbols<lenArray)
          rec_buffer[cntrFloatSymbols++] = *(float*)(data.data_ptr+k*sizeof(float));
    }
    std::cout<<"Symbols received: "<<(cntrFloatSymbols)<<'\n';
    if (cntrFloatSymbols==(lenArray))
    {
      std::cout<<"END OF RECIEVE DATA"<<std::endl;
      flagEndOfProc = true;
    }
    else
    {
      for (int cntrNop = 1;cntrNop<1000;cntrNop++)
        asm("nop");
    }

  }
  if (flagEndOfProc)
  { 
    strAnswer = "Recieve ";
    strAnswer.append(std::to_string(cntrFloatSymbols/2));
    strAnswer.append(" complex symbols. Mode = TX | Start TX (Freq = ");
    strAnswer.append(std::to_string(int(Freq*10)/10));
    strAnswer.append(" MHz. Fd = ");
    strAnswer.append(std::to_string(int(Fd)/1000));
    strAnswer.append(" MHz). Playing for ");
    strAnswer.append(std::to_string(int(SecondsToPlay)));
    strAnswer.append(" secodns. ");

  

    const char* tmp = strAnswer.c_str();
    bool statSend = client.sendData(tmp, sizeof(strAnswer)*strAnswer.length());
    statusOfSDR = playDeviceTx();
    strAnswer = "End streaming with status";
    strAnswer.append(std::to_string(statusOfSDR));
    const char* tmp2 = strAnswer.c_str();
    statSend = client.sendData(tmp2, sizeof(strAnswer)*strAnswer.length());
    if (statusOfSDR<0)
        strAnswer = "Please check you SDR. Server Ready To Next Command";
    else
        strAnswer = "Server Ready To Next Command";
    
    const char* tmp2_ = strAnswer.c_str();
    statSend = client.sendData(tmp2_, sizeof(strAnswer)*strAnswer.length());

/*
    char * arrayToTransfer = new char [cntrFloatSymbols*sizeof(float)];
    //client._status =  TcpServer::Client::status::up;
    statSend = client.sendData(arrayToTransfer, 1000);*/

  }
  cntr++;

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
