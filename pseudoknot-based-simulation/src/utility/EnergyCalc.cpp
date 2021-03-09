
#include "EnergyCalc.hpp"
#include <cassert>

namespace EnergyCalculator {

    double stemPairEnergy( const char i,const char j, const char p, const char q){ 
        bool NOT_CORRECT_PAIR=false;
    /*
    i p
    | |
    j q
    */
   int a,b;
   if(i == 'C' && j == 'G') a = 0;
   else if(i == 'G' && j == 'C') a = 1;
   else if(i == 'G' && j == 'U') a = 2;
   else if(i == 'U' && j == 'G') a = 3;
   else if(i == 'A' && j == 'U') a = 4;
   else if(i == 'U' && j == 'A') a = 5;
   else assert(NOT_CORRECT_PAIR);

   if(p == 'C' && q == 'G') b = 0;
   else if(p == 'G' && q == 'C') b = 1;
   else if(p == 'G' && q == 'U') b = 2;
   else if(p == 'U' && q == 'G') b = 3;
   else if(p == 'A' && q == 'U') b = 4;
   else if(p == 'U' && q == 'A') b = 5;
   else assert(NOT_CORRECT_PAIR);

    return EnergyConstants::stem[a][b];
    
    }


    double multiloopEnergy(int unpaired, int pairs) {
        double energy = EnergyConstants::multiloop[0] + EnergyConstants::multiloop[1]*unpaired + EnergyConstants::multiloop[2]*pairs;
        return energy;
    }


    double endPenalty(const char i, const char j){
        double energy = 0;
        bool NOT_CORRECT_PAIR=false;
        if (i == 'U'){
            if(j == 'A') energy += EnergyConstants::AU_end;
            else if (j == 'G') energy += EnergyConstants::GU_end;
            else assert(NOT_CORRECT_PAIR);
        } 
        else if (j == 'U'){
            if (i == 'A') energy += EnergyConstants::AU_end;
            else if (i == 'G') energy += EnergyConstants::GU_end; 
            else assert(NOT_CORRECT_PAIR);          
        }

        return energy;
            
    }

    double terminalMismatchEnergy(const char i, const char j, const char x, const char y) {

        bool NOT_CORRECT_PAIR=false;
        int a,c;
        if(i=='A' && j== 'U') a=0;
        else if(i=='C' && j== 'G') a=1;
        else if(i=='G' && j== 'C') a=2;
        else if(i=='G' && j== 'U') a=3;
        else if(i=='U' && j== 'A') a=4;
        else if(i=='U' && j== 'G') a=5;
        else assert(NOT_CORRECT_PAIR); 

        if(x=='A' && y=='A') c=0;
        else if(x=='A' && y=='C') c=1;
        else if(x=='A' && y=='G') c=2;
        else if(x=='A' && y=='U') c=3;
        else if(x=='C' && y=='A') c=4;
        else if(x=='C' && y=='C') c=5;
        else if(x=='C' && y=='G') c=6;
        else if(x=='C' && y=='U') c=7;
        else if(x=='G' && y=='A') c=8;
        else if(x=='G' && y=='C') c=9;
        else if(x=='G' && y=='G') c=10;
        else if(x=='G' && y=='U') c=11;
        else if(x=='U' && y=='A') c=12;
        else if(x=='U' && y=='C') c=13;
        else if(x=='U' && y=='G') c=14;
        else if(x=='U' && y=='U') c=15;
        else assert(NOT_CORRECT_PAIR); 

        return EnergyConstants::terminal_mismatch[a][c];

    }

// NOT USED (turner 2004)
    double hairpinSpecialEnergy(std::string str){
        

    if (str=="CAACG") return EnergyConstants::hairpinSpecial04[0];		   
    else if (str == "GUUAC") return EnergyConstants::hairpinSpecial04[1];	  
    else if (str == "CUACGG") return EnergyConstants::hairpinSpecial04[2];	  
    else if (str == "CUCCGG") return EnergyConstants::hairpinSpecial04[3];	  
    else if (str == "CUUCGG") return EnergyConstants::hairpinSpecial04[4];	  
    else if (str == "CUUUGG") return EnergyConstants::hairpinSpecial04[5];   
    else if (str == "CCAAGG") return EnergyConstants::hairpinSpecial04[6];	  
    else if (str == "CCCAGG") return EnergyConstants::hairpinSpecial04[7];	  
    else if (str == "CCGAGG") return EnergyConstants::hairpinSpecial04[8];	  
    else if (str == "CCUAGG") return EnergyConstants::hairpinSpecial04[9];	  
    else if (str == "CCACGG") return EnergyConstants::hairpinSpecial04[10];	  
    else if (str == "CCGCGG") return EnergyConstants::hairpinSpecial04[11];	  
    else if (str == "CCUCGG") return EnergyConstants::hairpinSpecial04[12];	  
    else if (str == "CUAAGG") return EnergyConstants::hairpinSpecial04[13];	  
    else if (str == "CUCAGG") return EnergyConstants::hairpinSpecial04[14];	  
    else if (str == "CUGCGG") return EnergyConstants::hairpinSpecial04[15];	  
    else if (str == "CAACGG") return EnergyConstants::hairpinSpecial04[16];	  
    else if (str == "ACAGUGCU") return EnergyConstants::hairpinSpecial04[17];	
    else if (str == "ACAGUGAU") return EnergyConstants::hairpinSpecial04[18];	
    else if (str == "ACAGUGUU") return EnergyConstants::hairpinSpecial04[19];	
    else if (str == "ACAGUACU") return EnergyConstants::hairpinSpecial04[20];	
    else return 0.0;
    }


    double internalLoop11Energy(const char i, const char j, const char p, const char q, const char x, const char y){
        int a,b,c;
        bool NOT_CORRECT_PAIR=false;
        //internal_loop11[a*6+b][c]
        if(i=='A' && j== 'U') a=0;
        else if(i=='C' && j== 'G') a=1;
        else if(i=='G' && j== 'C') a=2;
        else if(i=='U' && j== 'A') a=3;
        else if(i=='G' && j== 'U') a=4;
        else if(i=='U' && j== 'G') a=5;
        else assert(NOT_CORRECT_PAIR); 

        if(p=='A' && q== 'U') b=0;
        else if(p=='C' && q== 'G') b=1;
        else if(p=='G' && q== 'C') b=2;
        else if(p=='U' && q== 'A') b=3;
        else if(p=='G' && q== 'U') b=4;
        else if(p=='U' && q== 'G') b=5;
        else assert(NOT_CORRECT_PAIR); 

        if(x=='A' && y=='A') c=0;
        else if(x=='A' && y=='C') c=1;
        else if(x=='A' && y=='G') c=2;
        else if(x=='A' && y=='U') c=3;
        else if(x=='C' && y=='A') c=4;
        else if(x=='C' && y=='C') c=5;
        else if(x=='C' && y=='G') c=6;
        else if(x=='C' && y=='U') c=7;
        else if(x=='G' && y=='A') c=8;
        else if(x=='G' && y=='C') c=9;
        else if(x=='G' && y=='G') c=10;
        else if(x=='G' && y=='U') c=11;
        else if(x=='U' && y=='A') c=12;
        else if(x=='U' && y=='C') c=13;
        else if(x=='U' && y=='G') c=14;
        else if(x=='U' && y=='U') c=15;
        else assert(NOT_CORRECT_PAIR); 

        return EnergyConstants::internal_loop_11[a*6+b][c];
    }

    double danglingEnergy(const char i, const char j, const char x, bool prime3){
        int a,b;
        if(x=='A') b=0;
        else if (x=='C') b=1;
        else if (x=='G') b=2;
        else if (x=='U') b=3;

        if(i=='A' && j== 'U') a=0;
        else if(i=='C' && j== 'G') a=1;
        else if(i=='G' && j== 'C') a=2;
        else if(i=='G' && j== 'U') a=3;
        else if(i=='U' && j== 'A') a=4;
        else if(i=='U' && j== 'G') a=5;
        

        if(prime3) return EnergyConstants::dangling_3[a][b];
        else return EnergyConstants::dangling_5[a][b];

    }

    double internalmmEnergy(const char i, const char j, const char p, const char q){
        double energy =0;
    //GG is only considered in turner 2004

        if(i=='A' && j =='G')energy+= EnergyConstants::internalmm[0];        
        else if(i=='G' && j =='A')energy+= EnergyConstants::internalmm[1];  
    //    else if(i=='G' && j =='G')energy+= EnergyConstants::internalmm[2];
        else if(i=='U' && j =='U')energy+= EnergyConstants::internalmm[3];
        
        if(p=='A' && q =='G')energy+= EnergyConstants::internalmm[0];
        else if(p=='G' && q =='A')energy+= EnergyConstants::internalmm[1];
    //    else if(p=='G' && q =='G')energy+= EnergyConstants::internalmm[2];
        else if(p=='U' && q =='U')energy+= EnergyConstants::internalmm[3];

        return energy;
    }

    double hairpinTetraEnergy(std::string str){
        double energy =0;
        auto loop =str.substr(1,4);

        if(str[0] =='C'){
            assert(str.back() =='G');
            
            if(loop =="GGAA") energy= -3.0;
            else if(loop =="UUCG") energy= -3.0;
            else if(loop =="GUGA") energy= -3.0;
            else if(loop =="GAAA") energy= -3.0;
            else if(loop =="GCAA") energy= -3.0;
            else if(loop =="GAAG") energy= -2.5;
            else if(loop =="UACG") energy= -2.5;
            else if(loop =="GCGA") energy= -2.5;
            else if(loop =="GUAA") energy= -2.0;
            else if(loop =="UAAC") energy= -2.0;
            else if(loop =="GAGA") energy= -2.0;
            else if(loop =="GGGA") energy= -1.5;

        } else if (str[0] =='G'){
            assert(str.back() =='C' || str.back() == 'U');

            if(loop =="GGGA") energy= -3.0;
            else if(loop =="GUGA") energy= -3.0;
            else if(loop =="GAGA") energy= -3.0;
            else if(loop =="GAAA") energy= -3.0;
            else if(loop =="GCAA") energy= -2.5;
            else if(loop =="GAAG") energy= -1.5;
            else if(loop =="GGAG") energy= -1.5;
            else if(loop =="GCGA") energy= -1.5;
            else if(loop =="UGAA") energy= -1.5;
            else if(loop =="GGAA") energy= -1.5;

        } else if (str[0] =='U'){
            assert(str.back() =='A' ||str.back() =='G');

            if(str.back() =='G'){ 

                if(loop =="GAGA") energy= -2.5;
                else if(loop =="GAAA") energy= -2.0;

            }else{

                if(loop =="GAAA") energy= -1.5;
                else if(loop =="GGAA") energy= -1.5;
            }   

        } else if ( str[0]== 'A'){
            assert(str.back() == 'U');

            if(loop =="GAAA") energy= -2.0;
            else if(loop =="GCAA") energy= -1.5;
            else if(loop =="GUAA") energy= -1.5;
            else if(loop =="GUGA") energy= -1.5;
        }

        return energy;

    }

                                                

}; 