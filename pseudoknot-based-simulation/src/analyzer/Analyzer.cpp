
#include "Analyzer.hpp"
#include <fstream>
#include <omp.h>
#include <cstdio>
#include <cmath>
#include <vector>
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <functional>
#include <cassert>

void Analyzer::countAnalysis(    double starttme, double interval,double maxtime , int thr,std::vector<std::string> dotBrackets,std::string ifile,std::string ofile){
	std::cout<<"count Analysis " <<std::endl;
	 

		int datapoints=std::ceil((maxtime-starttme)/interval) +1;
		std::vector<int> all(datapoints*dotBrackets.size());
		std::string seq;
		int runsa=0;
		
		omp_set_num_threads(thr);
		
		#pragma omp parallel
		{
		std::ifstream infile2;
		std::string filenamein2 =  ifile+ std::to_string(omp_get_thread_num())+ ".txt";
		std::string line;
		std::string line2;
		
		int nst=0;
		infile2.open (filenamein2);
		
		
		std::vector<int> temp(datapoints*dotBrackets.size());
		int linenbmr=1;
		int runs=1;
		
		if( std::getline(infile2,line) && line[0]=='N'){
			if( omp_get_thread_num()==0){
				seq = line.substr(9);
			}
			linenbmr++;
			double tim=0;
			while(std::getline(infile2,line)){
				linenbmr++;
				if(line.find("New run")!=std::string::npos){
					linenbmr++;
					runs++;
					nst=0;
					if(runs%100==0 && omp_get_thread_num()==0)
					std::cout<<runs <<'\r'<< std::flush;
					if(!std::getline(infile2,line) ) break;
				}
				size_t pos = line.find("time: ");
				
				if(pos == std::string::npos){
					std::cout<< "error et line: " << linenbmr<< " no 'time: ' found" << std::endl;
					break; 
				}else{
					pos+=6;
									try{
						tim = std::stod(line.substr(pos));
					}catch(...){
						std::cout<< "error  line: " << linenbmr<< "thr: "<< omp_get_thread_num() << std::endl;
						break;
					}
					

				}


				if( tim>nst*interval+starttme &&nst*interval+ starttme<=maxtime ){
					for(int idx=0; idx< dotBrackets.size();idx++){
						// if(dotBrackets[idx]==line.substr(0,dotBrackets[idx].size())) temp[nst +idx*datapoints]++;
						if( line2.find(dotBrackets[idx]) != std::string::npos ) temp[nst +idx*datapoints]++;
					}
					nst++;
				}
				if(std::getline(infile2,line2))
				linenbmr++;
				else break;
			}
		}else std::cout<<"start error" <<std::endl;
		infile2.close();


		#pragma omp critical
		{
			runsa+=runs;
			for(int idx=0;idx<datapoints*dotBrackets.size(); idx++){
				all[idx] +=temp[idx];
			}
		}

		}
		

		std::cout<<"saving to file" <<std::endl;
		std::ofstream resultfile2;
		std::string filename2 =ofile+ ".txt";
		resultfile2.open (filename2);
		resultfile2 << seq << "\n";
        resultfile2 << "StartTime\tInterval\tMaxTime\n";
        resultfile2 << starttme <<"\t"<< interval <<"\t"<<maxtime <<"\n";
		resultfile2 << "runs\t" << runsa <<'\n';
		for(auto str : dotBrackets){
			resultfile2 << str << "\t";
		}
		resultfile2 << "\n";

		for(int idx=0;idx<datapoints; idx++){
			
			for(int idx2 =0;idx2<dotBrackets.size(); idx2++){
				resultfile2 << all[idx+ idx2*datapoints] << "\t";
			}
			resultfile2 << "\n";
		}

		resultfile2.close();
}
void Analyzer::countAnalysis(    double starttme, double interval,double maxtime , int thr,int mostFreq, std::string ifile,std::string ofile){
	assert(mostFreq>0);
	std::cout<<"count Analysis of "<< mostFreq << " most frequent structures" <<std::endl;
	 

		int datapoints=std::ceil((maxtime-starttme)/interval) +1;
		std::unordered_map<std::string,std::vector<int>> freqa;
		std::unordered_map<std::string,std::unordered_set<std::string>> familiesAll;
		
		int runsa=0;
		std::string seq;
		omp_set_num_threads(thr);
		
		#pragma omp parallel
		{
		std::ifstream infile2;
		std::string filenamein2 =  ifile+ std::to_string(omp_get_thread_num())+ ".txt";
		std::string line;
		std::string line2;
		
		int nst=0;
		infile2.open (filenamein2);
		//  int MFEnon=0; int MFEpk=0; int non1=0; int pk1=0;
		
		std::unordered_map<std::string,std::vector<int>> freqs;
		std::unordered_map<std::string,std::unordered_set<std::string>> families;
		int linenbmr=1;
		int runs=1;
		
		if( std::getline(infile2,line) && line[0]=='N'){
			if( omp_get_thread_num()==0){
				seq = line.substr(9);
			}
			linenbmr++;
			double tim=0;
			while(std::getline(infile2,line)){
				linenbmr++;
				if(line.find("New run")!=std::string::npos){
					linenbmr++;
					runs++;
					nst=0;
					if(runs%100==0 && omp_get_thread_num()==0)
					std::cout<<runs <<'\r'<< std::flush;
					if(!std::getline(infile2,line) ) break;
				}
				size_t pos = line.find("time: ");
				
				if(pos == std::string::npos){
					std::cout<< "error et line: " << linenbmr<< " no 'time: ' found" << std::endl;
					break; 
				}else{
					pos+=6;
									try{
						tim = std::stod(line.substr(pos));
					}catch(...){
						std::cout<< "error  line: " << linenbmr<< "thr: "<< omp_get_thread_num() << std::endl;
						break;
					}
					

				}


				if( tim>nst*interval+starttme &&nst*interval+ starttme<=maxtime ){
					 std::string line3 =line2;
					// line2 = line2.substr(26);
					if(freqs.count(line2)){
						freqs[line2][nst]+=1; 
						families[line2].insert(line3);
					}else{
						std::vector<int> vect(datapoints);
						vect[nst] =1;
						freqs[line2] = vect;
						families[line2]= std::unordered_set<std::string> {line3};
					}
					nst++;
				}
				if(std::getline(infile2,line2))
				linenbmr++;
				else break;
			}
		}else std::cout<<"start error" <<std::endl;
		infile2.close();


		#pragma omp critical
		{
			runsa+=runs;
			for( auto p: freqs){ 
				if(freqa.count(p.first)==1){
					for(int idx=0;idx<datapoints;idx++){
							freqa[p.first][idx]+=freqs[p.first][idx];
					}
					familiesAll[p.first].merge(families[p.first]);
						
				}else{
					
					freqa[p.first] = p.second;
					familiesAll[p.first] =families[p.first];
					
				}
			}
		}

		}
std::cout<< "different SCtrees "<<freqa.size()<< std::endl;
		auto cmp = [](std::pair<std::string,std::vector<int>> l,std::pair<std::string,std::vector<int>>r){ 
			auto v1= std::get<1>(l);
			auto v2= std::get<1>(r);
			int val1 =0;
			int val2 =0;
			for( int val :v1){
				val1+=val;
			}
			for( int val :v2){
				val2+=val;
			}
			return val1 > val2;
		};
		std::priority_queue<std::pair<std::string,std::vector<int>>,std::vector<std::pair<std::string,std::vector<int>> >, decltype(cmp) > que(cmp);

		for( auto e : freqa){
			que.push(e);
			if(que.size() >mostFreq) que.pop();
		}

		std::cout<<"saving to file" <<std::endl;
		std::ofstream resultfile2;
		std::string filename2 =ofile+ ".txt";
		resultfile2.open (filename2);
		resultfile2 << seq << "\n";
        resultfile2 << "StartTime\tInterval\tMaxTime\n";
        resultfile2 << starttme <<"\t"<< interval <<"\t"<<maxtime <<"\n";
		resultfile2 << "runs\t" << runsa <<'\n';
		std::vector<std::vector<int>> outvec;
		while(!que.empty()){
			resultfile2 << std::get<0>(que.top())<<" ("<<familiesAll[std::get<0>(que.top())].size() <<")\t";
			outvec.push_back(std::get<1>(que.top()));
			que.pop();
		}
		resultfile2 << '\n';
		int mf= std::min(mostFreq,(int)freqa.size());
		for(int idx=0;idx<datapoints; idx++){
			for(int idx2=0; idx2<mf-1;idx2++){
				resultfile2<< outvec[idx2][idx] <<"\t";
			}
			resultfile2<< outvec[mf-1][idx] <<"\n";
		}

		resultfile2.close();
}

void Analyzer::timeAnalysis(    double starttme, double interval,double maxtime , int thr,std::vector<std::string> dotBrackets, std::string ifile,std::string ofile){
	std::cout<<"time Analysis" <<std::endl;
	 

		int datapoints=std::ceil((maxtime-starttme)/interval) +1;
		std::vector<double> all(datapoints*dotBrackets.size());
		std::string seq;
		int runsa=0;
		
		omp_set_num_threads(thr);
		
		#pragma omp parallel
		{
		std::ifstream infile2;
		std::string filenamein2 =ifile + std::to_string(omp_get_thread_num())+ ".txt";
		std::string line;
		std::string line2;
		
		int nst=0;
		infile2.open (filenamein2);
		//  int MFEnon=0; int MFEpk=0; int non1=0; int pk1=0;
		
		std::vector<double> temp(datapoints*dotBrackets.size());
		int linenbmr=1;
		int runs=1;
		
		if( std::getline(infile2,line) && line[0]=='N'){
			if( omp_get_thread_num()==0){
				seq = line.substr(9);
			}
			linenbmr++;
			double tim=0;
			double timprev=0;
			while(std::getline(infile2,line)){
				linenbmr++;
				if(line.find("New run")!=std::string::npos){
					linenbmr++;
					runs++;
					nst=0;
					tim=0;
                    timprev=0;
					if(runs%100==0 && omp_get_thread_num()==0)
					std::cout<<runs << "\r" << std::flush;
					if(!std::getline(infile2,line) ) break;
				}
				size_t pos = line.find("time: ");
				
				if(pos == std::string::npos){
					std::cout<< "error et line: " << linenbmr<< " no 'time: ' found" << std::endl;
					break; 
				}else{
					pos+=6;
					try{
						tim = std::stod(line.substr(pos));
					}catch(...){
						std::cout<< "error  line: " << linenbmr<< "thr: "<< omp_get_thread_num() << std::endl;
						break;
					}
					

				}


				if(!line2.empty()){
                    if(tim>nst*interval+starttme  ){
							nst++;
                        }
                    if(nst*interval+ starttme<=maxtime){
                        

                        for(int idx=0; idx< dotBrackets.size();idx++){
							// if(dotBrackets[idx]==line2.substr(0,dotBrackets[idx].size())) temp[nst +idx*datapoints]+=( tim- timprev);
							if( line2.find(dotBrackets[idx]) != std::string::npos ) temp[nst +idx*datapoints]+=(tim- timprev);
						}
                    }
					 
					
					
				}
				if(std::getline(infile2,line2)){
					linenbmr++;
				}else break;

                timprev=tim;
				
			}
		}else std::cout<<"start error" <<std::endl;
		infile2.close();


		#pragma omp critical
		{
			runsa+=runs;
			for(int idx=0;idx<datapoints*dotBrackets.size(); idx++){
				all[idx] +=temp[idx];
			}
		}

		}


		std::cout<<"saving to file" <<std::endl;
		std::ofstream resultfile2;
		std::string filename2 = ofile + ".txt";
		resultfile2.open (filename2);
		resultfile2 << seq << "\n";
        resultfile2 << "StartTime\tInterval\tMaxTime\n";
        resultfile2 << starttme <<"\t"<< interval <<"\t"<<maxtime <<"\n";
		resultfile2 << "runs\t" << runsa <<'\n';
		for(auto str : dotBrackets){
			resultfile2 << str << "\t";
		}
		resultfile2 << "\n";

		for(int idx=0;idx<datapoints; idx++){
			
			for(int idx2 =0;idx2<dotBrackets.size(); idx2++){
				resultfile2 << all[idx+ idx2*datapoints] << "\t";
			}
			resultfile2 << "\n";
		}


		

		resultfile2.close();
}

std::pair<int,double>  Analyzer::neighboranalysis(RNASecondaryStructure sec){
	auto actions = sec.enumerateActions();
	double avE =0;
	for(auto a : actions){
		avE +=a.getEnergy();
	}

	return std::make_pair<int,double>(actions.size(),avE/actions.size());
}

