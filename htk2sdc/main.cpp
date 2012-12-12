#include <string>
#include <cmath>
#include <iostream>
#include <fstream>
#include <algorithm>
#include "tclap/CmdLine.h"

using namespace TCLAP;
using namespace std;

int _N;
int _d;
int _P;
int _k;
string _filelst;


int* pnSamples = new int;
int* psampPeriod = new int;
unsigned short* psampSize = new unsigned short;
unsigned short* pparmKind = new unsigned short;

float *HTKData;
float **t;

void convert_2(void* s);
void convert_4(void* s);
void convert_8(void* s);
void parseOptions(int argc, char** argv);
void dealingFileList(string filelist);
void readHTK(float* &HTKData,string filename);
void readHTKHeader(string filename);
void writeHTKheader(string filename,const int nSamples,const int sampPeriod,const unsigned short sampSize,const unsigned short parmKind);
void writeHTKdata(const float* htk,string filename,const int length);
void getSDCatFrame(const int N,const int d,const int P,const int k, float* &sdc,int frame);

int main(int argc, char** argv)
{
	parseOptions(argc,argv);
	cout<<"Dealing File:\""<<_filelst<<"\" with Parameter:\tN:"<<_N<<"    "<<"d:"<<_d<<"    "<<"P:"<<_P<<"    "<<"k:"<<_k<<endl;
    
    dealingFileList(_filelst);
}

void dealingFileList(string flst){
    ifstream fileList(flst.c_str(),ios::in);
    
    if(!fileList.is_open()){
        cout<<"\tERROR filelist does not exist:"<<flst<<endl;
        exit(EXIT_FAILURE);
    }
    
    string infile;
    string outfile;
    while(!fileList.eof()){
        fileList>>infile;
        fileList>>outfile;
        cout<<infile<<" -> "<<outfile<<"\t"<<flush;
        
        readHTKHeader(infile);
        
        HTKData = (float*) malloc((*pnSamples)*(*psampSize)*sizeof(float));
        readHTK(HTKData,infile);
        writeHTKheader(outfile, *pnSamples, *psampPeriod, _k*_N, *pparmKind);
        
        t = (float**)malloc(*psampSize*sizeof(float*));
        for(int i=0;i<*psampSize;i++)
            t[i] = (float*)malloc(*pnSamples*sizeof(float));
        
        //convert raw data to 2-d matrix
        for(int i=0;i<*pnSamples**psampSize;i++)
            t[i%*psampSize][i/ *psampSize] = HTKData[i];
        
        int old = 0;
        for(int i=0;i<*pnSamples;i++){
            int now = (i*100/ *pnSamples);
            float* SDC;
            getSDCatFrame(_N, _d, _P, _k, SDC, i);
            writeHTKdata(SDC, outfile, _k*_N);
            
            if(now-old==10){
                cout<<"."<<flush;
                old+=10;
            }
            delete SDC;
        }
        
        for(int i=0; i<*psampSize; i++)
            delete [] t[i];
        
        delete[] t;
        delete[] HTKData;
        cout<<"done!"<<endl;
    }
    fileList.close();
    cout<<"\n\nAll these jobs has been done!"<<endl;
    
}

void readHTKHeader(string filename){
    ifstream HTKfile(filename.c_str(),ios::binary);
    
    if(!HTKfile.is_open()){
        cout<<"\tERROR file does not exist:"<<filename<<endl;
        exit(EXIT_FAILURE);
    }
    
	HTKfile.read((char*)pnSamples,sizeof(*pnSamples));
	HTKfile.read((char*)psampPeriod,sizeof(*psampPeriod));
	HTKfile.read((char*)psampSize,sizeof(*psampSize));
	HTKfile.read((char*)pparmKind,sizeof(*pparmKind));
    
	convert_4(pnSamples);
	convert_4(psampPeriod);
	convert_2(psampSize);
	convert_2(pparmKind);
    
	*psampSize/=4;
}

void getSDCatFrame(const int N,const int d,const int P,const int k, float* &sdc,int frame){
    int maxk = (*pnSamples-1-2*d)/P+1;
    if(k>maxk){
        cout<<"\tERROR: parameter 'k' is overflow"<<endl;
        exit(EXIT_FAILURE);
    }
    if(N>*pnSamples){
        cout<<"\tERROR: parameter 'N' is overflow"<<endl;
        exit(EXIT_FAILURE);
    }
    if(frame>*pnSamples){
        cout<<"\tERROR: frame is overflow"<<endl;
        exit(EXIT_FAILURE);
    }

    
    float **t2;
    t2 = (float**)malloc(N*sizeof(float*));
    for(int i=0;i<N;i++)
        t2[i] = (float*)malloc(k*sizeof(float));
    
    sdc = (float*)malloc(N*maxk);
    
    for(int i=frame;(i-frame)/P<k;i+=P){
        float tmp= 0;
        for(int j=0;j<N;j++){
            if(i-d>=0&&i+d<*pnSamples){
                tmp = t[j][i+d]-t[j][i-d];
            }
            t2[j][(i-frame)/P] = tmp;
        }
    }
    
    //convert 1frame-sdc 2d-martix to 1d-martix
    for(int i=0;i<k;i++){
        for(int j=0;j<N;j++){
            float tmp = t2[j][i];
            sdc[i*N+j] = tmp;
        }
    }
    
    for(int i=0;i<N;i++)
        delete[] t2[i];
    delete[] t2;
}

void writeHTKheader(string filename,const int nSamples,const int sampPeriod,const unsigned short sampSize,const unsigned short parmKind){
    ofstream HTKfile(filename.c_str(),ios::binary);
   
    
    int tnSamples = nSamples;
    int tsampPeriod = sampPeriod;
    unsigned short tsampSize = sampSize;
    unsigned short tparamKind = parmKind;
    tsampSize*=4;
    convert_4(&tnSamples);
    convert_4(&tsampPeriod);
    convert_2(&tsampSize);
    convert_2(&tparamKind);
    
    HTKfile.write((char*)&tnSamples,sizeof(tnSamples));
    HTKfile.write((char*)&tsampPeriod,sizeof(tsampPeriod));
    HTKfile.write((char*)&tsampSize,sizeof(tsampSize));
    HTKfile.write((char*)&tparamKind,sizeof(tparamKind));
}

void writeHTKdata(const float* htk,string filename,const int length){
    ofstream HTKfile(filename.c_str(),ios::binary|ios::app);
    for(int i=0;i<length;i++){
        float tmp = htk[i];
        convert_4(&tmp);
        HTKfile.write((char*)&tmp,sizeof(tmp));
    }
}

void readHTK(float* &HTKData,string filename){
	ifstream HTKfile(filename.c_str(),ios::binary);
    
    if(!HTKfile.is_open()){
        cout<<"\tERROR file does not exist:"<<filename<<endl;
        exit(EXIT_FAILURE);
    }
    
	HTKfile.read((char*)pnSamples,sizeof(*pnSamples));
	HTKfile.read((char*)psampPeriod,sizeof(*psampPeriod));
	HTKfile.read((char*)psampSize,sizeof(*psampSize));
	HTKfile.read((char*)pparmKind,sizeof(*pparmKind));
    
	convert_4(pnSamples);
	convert_4(psampPeriod);
	convert_2(psampSize);
	convert_2(pparmKind);
    
	*psampSize/=4;
    
    for(int i=0;i<(*psampSize)*(*pnSamples);i++){
        HTKfile.read((char*)&HTKData[i],sizeof(HTKData[i]));
        convert_4(&HTKData[i]);
    }
}

void convert_2(void* s){
	unsigned short *x = (unsigned short*)s;
	*x = ((*x & 0xff) << 8) | (*x >> 8);
}

void convert_4(void* s){
	unsigned int *x = (unsigned int*)s;
	*x = (*x << 24) |
    ((*x << 8) & 0x00ff0000u) |
    ((*x >> 8) & 0x0000ff00u) |
    (*x >> 24);
}

void convert_8(void* s){
	unsigned long long int *x = (unsigned long long int *)s;
	*x = (*x<<56)|
    ((*x<<40)	&	0x00ff000000000000u)|
    ((*x<<24)	&	0x0000ff0000000000u)|
    ((*x<<8)	&	0x000000ff00000000u)|
    ((*x>>8)	&	0x00000000ff000000u)|
    ((*x>>24)	&	0x0000000000ff0000u)|
    ((*x>>40)	&	0x000000000000ff00u)|
    ((*x>>56));
}

void parseOptions(int argc, char** argv){
	try {
		CmdLine cmd("htk2sdc", ' ', "0.00001" );
        
		//
		// Define arguments
		//
        
		ValueArg<string> arg_filelist("f","fileList","Filelist to convert", true, "", "string");
		cmd.add(arg_filelist);
        
		ValueArg<int> arg_N("N","parameter-N","The number of cepstral coefficients computed at each frame.",true,0,"int");
		cmd.add(arg_N);
        
		ValueArg<int> arg_d("d","parameter-d","The time advance and delay for the delta computation.",true,0,"int");
		cmd.add(arg_d);
        
		ValueArg<int> arg_P("P","parameter-P","The time shift between consecutive blocks.",true,0,"int");
		cmd.add(arg_P);
        
		ValueArg<int> arg_k("k","parameter-k","The number of blocks whose delta coefficients are concatenated to form the final feature vector",true,0,"int");
		cmd.add(arg_k);
        
		//
		// Parse the command line.
		//
		cmd.parse(argc,argv);
        
		//
		// Set variables
		//
		_filelst = arg_filelist.getValue();
		_N = arg_N.getValue();
		_d = arg_d.getValue();
		_P = arg_P.getValue();
		_k = arg_k.getValue();
	} catch ( ArgException& e ){ 
		cout << "\tERROR: " << e.error() << " " << e.argId() << endl; 
	}
}
