#include "Simulator.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <random>
#include <algorithm>
#include <cassert>
#include "constants/Constants.hpp"
#include<omp.h>
#include <chrono> 

Simulator::Simulator() {}

void Simulator::runSim(double duration, RNASecondaryStructure secondary,std::string file ,bool cotranscription, double transSpeed){

    
    {
    double tautrans = 1.0/transSpeed;
    std::ofstream resultfile;
    std::string filename = file;
    filename += std::to_string(omp_get_thread_num());
    filename += ".txt";
    resultfile.open (filename ,std::ios::app );
    
    resultfile << "New run: "<< secondary.getSequence() <<"\n";
    
     
      
    std::mt19937 gen(static_cast<long unsigned int>(std::chrono::high_resolution_clock::now().time_since_epoch().count()) +omp_get_thread_num()); 
    
    double energy =secondary.evaluateEnergy();
    resultfile << "initial energy: " << energy << ", time: 0.0\n";
    resultfile << secondary.getDotBracket() <<"\n";
    double time = 0.0;
    double sumtrans = tautrans;
    bool trans = true;
    // std::cout<<"time: 0.0"<< std::endl;
    while(time < duration){
    
        std::vector<Action> act= secondary.enumerateActions();
        //################################## debugging
        // std::vector<Action> act2= secondary.enumerateActionsOrder();
        // std::sort(act.begin(),act.end());
        // std::sort(act2.begin(),act2.end());
        // assert(act.size()==act2.size());
        // for(int ii =0; ii<act.size(); ii++) assert(act[ii]==act2[ii]);
        //################################## debugging
        assert(act.size() >0 || cotranscription); //not to be asserted in cotranscriptional fold
        // counting cumulative rates
        double taufold=std::numeric_limits<double>::max();
        std::vector<double> actidx(act.size());
        std::exponential_distribution timeDis;
        if(!act.empty()){
            
            actidx[0]= std::exp(-act[0].getEnergy()/ (EnergyConstants::R * EnergyConstants::T *2) );
            for(int idx=1;idx<(int)act.size(); idx++){
                actidx[idx]=actidx[idx-1] + std::exp(-act[idx].getEnergy()/ (EnergyConstants::R * EnergyConstants::T *2) );
            }
            std::exponential_distribution timeDis(actidx.back()*EnergyConstants::k0);
            taufold= timeDis(gen);
        }
        
        


        // cotranscriptional folding
        
        if(cotranscription)
        while(sumtrans < taufold+time && trans){
            time = sumtrans;
            trans =secondary.increaseFoldingLength();
            if(!trans) break;
            auto addAct = secondary.additonalActions();
            act.insert(act.end(),addAct.begin(),addAct.end());
            if(!addAct.empty()){
                if(actidx.empty()){
                    actidx.push_back( std::exp(-addAct[0].getEnergy()/ (EnergyConstants::R * EnergyConstants::T *2) ));
                }
                for(int idx=1;idx<(int)addAct.size(); idx++ ){
                    actidx.push_back( actidx.back() + std::exp(-addAct[idx].getEnergy()/ (EnergyConstants::R * EnergyConstants::T *2) ));
                }
            }
            
            if(!act.empty()){
                timeDis = std::exponential_distribution(actidx.back()*EnergyConstants::k0);
                taufold = timeDis(gen);
            }
            
            sumtrans += tautrans;


            // std::cout<<"time: "<<time<< " Energy: "<< energy << " thr: "<< omp_get_thread_num()<< std::endl;
            // std::cout<<secondary.getDotBracket()<<"    thr: "<< omp_get_thread_num()<< std::endl;
            resultfile << "energy: " << energy << ", time: "<<time<< " fl: " << secondary.getFoldingLength() << "\n" ;
            
            resultfile << secondary.getDotBracket()<< "\n";
        }

        time += taufold;
        std::uniform_real_distribution<> dis(0.0, actidx.back());
        double find =dis(gen);
        auto actionidx =std::upper_bound(actidx.begin(),actidx.end(),find);
        int idx =actionidx-actidx.begin();
        secondary.executeAction(act[idx]);
        //double e = secondary.getEnergy();
        energy =secondary.evaluateEnergy(); // optimize to getEnergy()
        //assert(e==energy);
        // std::cout<<"time: "<<time<< " Energy: "<< energy << " thr: "<< omp_get_thread_num()<< std::endl;
        // std::cout<<secondary.getDotBracket()<<"    thr: "<< omp_get_thread_num()<< std::endl;
        resultfile << "energy: " << energy << ", time: "<<time/* << " action:" << ((act[idx].getType()==Addition) ? " Addition, ": " Deletion, ") ;
        resultfile <<"i: " <<act[idx].i<< " j: " <<act[idx].j<< " \u0394E " <<act[idx].getEnergy() */ << "\n";
        resultfile << secondary.getDotBracket() <<"\t" <<  secondary.getSCTree() <<"\n";
        // std::cout << secondary.getDotBracket() <<"\t" <<std::endl;
        // std::cout <<  secondary.getSCTree() << std::endl;

    }


    resultfile.close();
    }
    
}

void Simulator::simSteps(int maxSteps, RNASecondaryStructure secondary,std::string file ,bool cotranscription, double transSpeed){

    
    {
    double tautrans = 1.0/transSpeed;
    std::ofstream resultfile;
    std::string filename = file;
    filename += std::to_string(omp_get_thread_num());
    filename += ".txt";
    resultfile.open (filename ,std::ios::app );
    
    resultfile << "New run: "<< secondary.getSequence() <<"\n";
    
     
      
    std::mt19937 gen(static_cast<long unsigned int>(std::chrono::high_resolution_clock::now().time_since_epoch().count()) +omp_get_thread_num()); 
    
    double energy =secondary.evaluateEnergy();
    resultfile << "initial energy: " << energy << ", time: 0.0\n";
    resultfile << secondary.getDotBracket() <<"\n";
    double time = 0.0;
    double sumtrans = tautrans;
    bool trans = true;
    // std::cout<<"time: 0.0"<< std::endl;
    int steps=0;
    while(steps<maxSteps){
    
        std::vector<Action> act= secondary.enumerateActions();
        //################################## debugging
        // std::vector<Action> act2= secondary.enumerateActionsOrder();
        // std::sort(act.begin(),act.end());
        // std::sort(act2.begin(),act2.end());
        // assert(act.size()==act2.size());
        // for(int ii =0; ii<act.size(); ii++) assert(act[ii]==act2[ii]);
        //################################## debugging
        assert(act.size() >0 || cotranscription); //not to be asserted in cotranscriptional fold
        // counting cumulative rates
        double taufold=std::numeric_limits<double>::max();
        std::vector<double> actidx(act.size());
        std::exponential_distribution timeDis;
        if(!act.empty()){
            
            actidx[0]= std::exp(-act[0].getEnergy()/ (EnergyConstants::R * EnergyConstants::T *2) );
            for(int idx=1;idx<(int)act.size(); idx++){
                actidx[idx]=actidx[idx-1] + std::exp(-act[idx].getEnergy()/ (EnergyConstants::R * EnergyConstants::T *2) );
            }
            std::exponential_distribution timeDis(actidx.back()*EnergyConstants::k0);
            taufold= timeDis(gen);
        }
        
        


        // cotranscriptional folding
        
        if(cotranscription)
        while(sumtrans < taufold+time && trans){
            time = sumtrans;
            trans =secondary.increaseFoldingLength();
            if(!trans) break;
            auto addAct = secondary.additonalActions();
            act.insert(act.end(),addAct.begin(),addAct.end());
            if(!addAct.empty()){
                if(actidx.empty()){
                    actidx.push_back( std::exp(-addAct[0].getEnergy()/ (EnergyConstants::R * EnergyConstants::T *2) ));
                }
                for(int idx=1;idx<(int)addAct.size(); idx++ ){
                    actidx.push_back( actidx.back() + std::exp(-addAct[idx].getEnergy()/ (EnergyConstants::R * EnergyConstants::T *2) ));
                }
            }
            
            if(!act.empty()){
                timeDis = std::exponential_distribution(actidx.back()*EnergyConstants::k0);
                taufold = timeDis(gen);
            }
            
            sumtrans += tautrans;


            // std::cout<<"time: "<<time<< " Energy: "<< energy << " thr: "<< omp_get_thread_num()<< std::endl;
            // std::cout<<secondary.getDotBracket()<<"    thr: "<< omp_get_thread_num()<< std::endl;
            resultfile << "energy: " << energy << ", time: "<<time<< " fl: " << secondary.getFoldingLength() << "\n" ;
            
            resultfile << secondary.getDotBracket()<< "\n";
        }

        time += taufold;
        steps++;
        std::uniform_real_distribution<> dis(0.0, actidx.back());
        double find =dis(gen);
        auto actionidx =std::upper_bound(actidx.begin(),actidx.end(),find);
        int idx =actionidx-actidx.begin();
        secondary.executeAction(act[idx]);
        //double e = secondary.getEnergy();
        energy =secondary.evaluateEnergy(); // optimize to getEnergy()
        //assert(e==energy);
        // std::cout<<"time: "<<time<< " Energy: "<< energy << " thr: "<< omp_get_thread_num()<< std::endl;
        // std::cout<<secondary.getDotBracket()<<"    thr: "<< omp_get_thread_num()<< std::endl;
        resultfile << "energy: " << energy << ", time: "<<time/* << " action:" << ((act[idx].getType()==Addition) ? " Addition, ": " Deletion, ") ;
        resultfile <<"i: " <<act[idx].i<< " j: " <<act[idx].j<< " \u0394E " <<act[idx].getEnergy() */ << "\n";
        resultfile << secondary.getDotBracket() <<"\t" <<  secondary.getSCTree() <<"\n";
        // std::cout << secondary.getDotBracket() <<"\t" <<std::endl;
        // std::cout <<  secondary.getSCTree() << std::endl;

    }


    resultfile.close();
    }
    
}


