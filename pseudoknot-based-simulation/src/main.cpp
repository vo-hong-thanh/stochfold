//the main file 

#include "primary/RNAPrimarySequence.hpp"
#include "secondary/RNASecondaryStructure.hpp"
#include <iostream>
#include "simulator/Simulator.hpp"
#include <fstream>
#include <omp.h>
#include <cstdio>
#include <cmath>
#include "analyzer/Analyzer.hpp"
#include <random>
#include <sys/time.h>
#include <ctime>
#include <unistd.h>

static double get_time() {
    struct timeval tm;
    gettimeofday(&tm, NULL);
    return static_cast<double>(tm.tv_sec)
        + static_cast<double>(tm.tv_usec) / 1E6;
}


void printExample(){
	std::ifstream infile;
	std::string filenamein ="../exampleinput.txt";
	std::string line;
	infile.open (filenamein);
	while(std::getline(infile,line)){
		std::cout<<line<<std::endl;
	}
	infile.close();
}



std::string gen_random(const int len) {
    
    std::string tmp_s;
    static const char alphanum[] =
        "ACGU";
    
    srand( (unsigned) time(NULL) * getpid());

    tmp_s.reserve(len);

    for (int i = 0; i < len; ++i) 
        tmp_s += alphanum[rand() % (sizeof(alphanum) - 1)];
    
    
    return tmp_s;
    
}



int main() {
	

 

	
	
	//pseudoknot numbering works for ~all seq
	std::vector<int> fig2pk{0,0,1,1,2,2,3,3,4,4};

   std::vector<std::string> dbss{
".((((..((((.((((((....))).))))).))..).)))......",
"..((((((((......))))))))....(((.(....).))).....",
"..((((((((......))))))))..(((((.........)).))).",
"..((((((((......))))))))....(((........))).....",
"..((((((((......))))))))...(((([[[[....))))]]]]",
".(((((...(.(((.((.((....)))).))).)..)))))......",
"..((((((((......))))))))...((((((...)).))))....",
"..(((((((([[[[[.))))))))...((((.]]]]]..))))....",
"...(((((((......)))))))....((((.(....).))))....",
"...(((((((......)))))))....((((........))))....",
"..((((((((......))))))))...((((.(....).))))....",
"..((((((((......))))))))...((((........))))....",
"(((.((((((.[[[[.)))))).....((((.]]]]...)))))))."
   };
	
	// for(int idx=0;idx<13;idx++){
int idx=1;
	
	RNAPrimarySequence p = RNAPrimarySequence( "CGGUAGCGCGAACCGUUAUCGCGCA" /* "GCCGGGCGCGGUGGCGCGUGCCUGUAGUCCCAGCUACUCGGGAGGCU"*/);
	RNASecondaryStructure sec = RNASecondaryStructure(p, "(....[[[[[[...)...]]]]]]." /* dbss[idx]*/,fig2pk);
   
	
    // std::cout << sec.getSCTree() << std::endl;
	// std::cout << sec.evaluateEnergy() << std::endl;
	// }

	std::string fileid = "FILEID";
    // thridx +".txt" is added to filenames
    std::string sofile="results/raw/results" + fileid;
    std::string aofile="results/distributions/disCounts"+fileid;
    // std::string ofile="results/distributions/disTime"+fileid;
	int thr=omp_get_max_threads();
	int runs=0;
	double duration = 0;
	double start=0;
	double end=30;
	double interval=1;
	int mostfreq=8;
	// std::vector<std::string> dbs= {".....(((((((....).))))))."	,"......(((((.......)))))..",".(((((((.....))))))).....",".(((.[[[[[[)))....]]]]]]."	,	"((((((((.....))))))))....",	"((((.[[[[[.))))....]]]]]."	,	"(((..[[[[[[.)))...]]]]]]."	,	".....(((((.........)))))."};
	//  std::vector<std::string> dbs= {"((((.[[[[[[))))...]]]]]]."};
	// std::vector<std::string> conformations = {".....((((((.......)))))).","(....[[[[[[...)...]]]]]].","...(.[[[[[[)......]]]]]].","((((.[[[[[[))))...]]]]]].", "......(((((.......))))).."	,	".((((((.......))))))....."	,	".(((((((.....)))))))....."	,	".(((.[[[[[[)))....]]]]]]."	,	"((((((((.....))))))))...."	,	"((((.[[[[[.))))....]]]]]."	,	"(((..[[[[[[.)))...]]]]]].",".....(((((.........)))))."};
	std::vector<std::string> conformations;










	
	std::ifstream infile;
	std::string filenamein ="input.txt";
	std::string line;
	infile.open (filenamein);


	if(!std::getline(infile,line) || line[0] != '>'){
	    std::cout<< "NOT VALID INPUT, first character shlould be '>'" << std::endl;
		
	    return 0;
	}else if(std::getline(infile,line) && line.size() >0){
	    line.erase(std::remove(line.begin(), line.end(), ' '), line.end());
	    for( auto c : line){
		    if( !(c== 'A' || c== 'U' || c== 'G' || c== 'C') ){
			    std::cout<< "NOT VALID character: " <<c  << std::endl;
			    return 0;
		    }
	    }
	    p= RNAPrimarySequence(line, 1 ); // <- add initial folding Len here
	}else {
		printExample();
		return 0;
	}
	if(std::getline(infile,line)){
	    sec = RNASecondaryStructure(p,line,fig2pk);
	    std::cout<< line  << std::endl;
	}else {
		printExample();
		return 0;
	}
	std::vector<bool> params(12);
	// std::vector<std::string> paramwords{"fileid=" ,"sofile=","aofile=", "thr=", "runs=","mostfreq=","duration=","start=","interval=","end=","dotBrackets={", "timeAnalysis"  }
	int pos=0;
	int pos1=-1;
	while(std::getline(infile,line) && line != "---"){
		


		
		if(line[0] != '#'){
			if(!params[0]){
				pos =line.find("fileid=");
				if(pos!= std::string::npos){
					pos1=line.find(",",pos);
					fileid=line.substr(pos+7,pos1-pos);
					params[0] =true;
				}
			}
			if(!params[1]){
				pos =line.find("sofile=");
				if(pos!= std::string::npos){
					pos1=line.find(",",pos);
					sofile=line.substr(pos+7,pos1-pos);
					params[1] =true;
				}
			}
			if(!params[2]){
				pos =line.find("aofile=");
				if(pos!= std::string::npos){
					pos1=line.find(",",pos);
					aofile=line.substr(pos+7,pos1-pos);
					params[2] =true;
				}
			}
			if(!params[3]){
				pos =line.find("thr=");
				if(pos!= std::string::npos){
					pos1=line.find(",",pos);
					try{
						thr= std::stoi( line.substr(pos+4,pos1-pos));
						params[3] =true;
					}catch(...){ }
					
				}
			}
			if(!params[4]){
				pos =line.find("runs=");
				if(pos!= std::string::npos){
					pos1=line.find(",",pos);
					try{
						runs= std::stoi( line.substr(pos+5,pos1-pos));
						params[4] =true;
					}catch(...){ }
				}
					
			}
			if(!params[5]){
				pos =line.find("mostfreq=");
				if(pos!= std::string::npos){
					pos1=line.find(",",pos);
					try{
						mostfreq= std::stoi( line.substr(pos+9,pos1-pos));
						params[5] =true;
					}catch(...){ }
					
				}
			}
			if(!params[6]){
				pos =line.find("duration=");
				if(pos!= std::string::npos){
					pos1=line.find(",",pos);
					try{
						duration= std::stod( line.substr(pos+9,pos1-pos));
						params[6] =true;
					}catch(...){ }
					
				}
			}
			if(!params[7]){
				pos =line.find("start=");
				if(pos!= std::string::npos){
					pos1=line.find(",",pos);
					try{
						start= std::stod( line.substr(pos+6,pos1-pos));
						params[7] =true;
					}catch(...){ }
					
				}
			}
			if(!params[8]){
				pos =line.find("interval=");
				if(pos!= std::string::npos){
					pos1=line.find(",",pos);
					try{
						interval= std::stod( line.substr(pos+9,pos1-pos));
						params[8] =true;
					}catch(...){ }
					
				}
			}
			if(!params[9]){
				pos =line.find("end=");
				if(pos!= std::string::npos){
					pos1=line.find(",",pos);
					try{
						end= std::stod( line.substr(pos+4,pos1-pos));
						params[9] =true;
					}catch(...){ }
					
				}
			}
			if(!params[10]){
				pos =line.find("conformations={");
				if(pos!= std::string::npos){
					while(std::getline(infile,line) && line != "}"){
						if(line[0]!='#')
						conformations.push_back(line);
					}
					if (!conformations.empty()) params[10]=true;
				}
			}
			if(!params[11]){
				pos =line.find("timeAnalysis");
				if(pos!= std::string::npos){
					params[11] =true;
				}
			}

			
		}
	}
	infile.close();

	sofile = sofile + fileid;
	aofile = aofile + fileid;
	



	
	omp_set_num_threads(thr);
	

    
	///timetesting
	// auto lens = std::vector<int>{75,200,300};
	// runs=1;
	// for(int len : lens)
	// for(int dummy=0;dummy<10;dummy++){
	// 	RNAPrimarySequence p1 = RNAPrimarySequence(gen_random(len));
	// 	RNASecondaryStructure sec = RNASecondaryStructure(p1,std::string(len,'.'),fig2pk);

	// //timetesting
	double t0 = get_time();
	if(params[1] && params[4] && params[6]) //uncomment after timetest
	#pragma omp parallel
	{
	
	
	  Simulator sim;
	  for(int idx=0;idx<runs;idx++){
        if(omp_get_thread_num() ==3)
	    std::cout<< "  run: " << idx +1<< "/"<< runs << '\r' <<std::flush;
        //for cotranscriptional add "true" as 4th parameter and add folding length to primary
		sim.runSim(duration,sec,sofile ,true ,200.0); //uncomment timetest
		// sim.simSteps(1000,sec,sofile ,false ); //timetest
	  }
	
	}
	double t1 = get_time();
	// std::cout << "len: " << len << "iter: " << dummy; // timetest
	std::cout << " sim_took " << t1 - t0 << " seconds" << std::endl;
	// } //timetesting
	// return 0; // timetesting
	

	Analyzer ana;
	
	if(params[1] && params[2] && params[7] &&params[8] && params[9]){

		if(params[5])
		ana.countAnalysis(start,interval,end,thr,mostfreq,sofile,aofile);
		else if(params[10]){
			if(params[11])
			ana.timeAnalysis(start,interval,end,thr,conformations,sofile,aofile);
			else
			ana.countAnalysis(start,interval,end,thr,conformations,sofile,aofile);
		}else {
		printExample();
		return 0;
	}
		
	}else {
		printExample();
		return 0;
	}
	

	// auto pa =ana.neighboranalysis(sec);
	// std::cout << std::get<0>(pa) <<" "<< std::get<1>(pa) << std::endl;
	std::cout << "END OF TEST" << std::endl;
	
}

