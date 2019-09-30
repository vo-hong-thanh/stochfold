/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package utility;

import java.util.Random;
import model.RNAPrimarySequence;
import constants.ENERGY_CONSTANTS;
import constants.NUCLEOTIDE_CODE;
import mfe_optimization.ZukerAlgorithm;

//1 Kilocalories = 4184 Jules
/**
 *
 * @author vot2
 */
public class ComputingMachine {
    //compute firing time
    public static double computeTime(Random rand, double totalRate) {
        double delta = Double.MAX_VALUE;
        if (totalRate != 0.0) {
            delta = Math.log(1 / rand.nextDouble()) / totalRate;
        }
        return delta;
    }

    //compute rate
    public static double computeRate(double energy) {
        return Math.exp(-energy / (2 * ENERGY_CONSTANTS.GAS_CONSTANT * ENERGY_CONSTANTS.T));
    }

    //eveluate optimal energy
    public static double computeOptimalEnergy(RNAPrimarySequence sequence) {
        ZukerAlgorithm zuker = new ZukerAlgorithm();
        zuker.config(sequence);
        zuker.evaluateOptimalEnergy();
        return zuker.getOptimalEnergy();
    }
    
    //eveluate energy
    public static double computeStructureEnergy(RNAPrimarySequence sequence, int[] pairList, int currentFoldingLength) {
        double energy = 0.0;
        double ee;
        
        for (int i = 0; i < currentFoldingLength; i++) {
            if (pairList[i] > i) {
//                System.out.println("Consider component enclosed by the base-pair (" + i + ", " + pairList[i] + ")");
                ee = componentEnergy(sequence, pairList, i);
                energy += ee;

//                ee = exteriorloopEnergy(i, pairList[i], sequence, pairList, currentFoldingLength);
//                System.out.println(" Exterior componenet energy " + ee + " (kcal)");
//                energy += ee;
                i = pairList[i];
            }
        }
                
        return energy;
    }    
   
//    private static double exteriorloopEnergy(int fivePrime, int threePrime, RNAPrimarySequence sequence, int[] pairList, int currentFoldingLength) {
//        double energy = 0.0;
//        double ee;
//        
////        System.out.print(" Exterior loop (");
//        char iNucleotide = sequence.getNucleotide(fivePrime);
//        char jNucleotide = sequence.getNucleotide(threePrime);
//
////        if (checkAUEnd(iNucleotide, jNucleotide)) {
////            ee = ENERGY_CONSTANTS.AU_end;
//////                    System.out.print(" AU end penalty by ( i = " + fivePrime + " , j = " + threePrime + " ) = " + ee + " (kcal),");
////            energy += ee;
////        } else if (checkGUEnd(iNucleotide, jNucleotide)) {
////            ee = ENERGY_CONSTANTS.GU_end;
//////                    System.out.print(" GU end penalty by ( i = " + fivePrime + " , j = " + threePrime + " ) = " + ee + " (kcal),");
////            energy += ee;
////        }
//
//        //terminal mismatch
//        if (fivePrime > 0 && (threePrime + 1 < currentFoldingLength) &&  pairList[threePrime + 1] == -1) {
//            ee = lookupTerminalMismatchEnergy(iNucleotide, jNucleotide, sequence.getNucleotide(fivePrime - 1), sequence.getNucleotide(threePrime + 1));
////                    System.out.print(" Terminal mismatch = " + ee + " (kcal)");
//            energy += ee;
//        }
//        //5' dangling end        
//        else if (fivePrime > 0) {
//            ee = lookupDanglingEnergy(sequence.getNucleotide(fivePrime), sequence.getNucleotide(pairList[fivePrime]), sequence.getNucleotide(fivePrime - 1), true);
////                    System.out.println(" 5' dangling end = " + ee + " (kcal)");
//            energy += ee;
//        }//3' dangling end
//        else if ((threePrime + 1 < currentFoldingLength) &&  pairList[threePrime + 1] == -1) 
//        {
//            ee = lookupDanglingEnergy(sequence.getNucleotide(pairList[threePrime]), sequence.getNucleotide(threePrime), sequence.getNucleotide(threePrime + 1), false);
////                    System.out.println(" 3' dangling end = " + ee + " (kcal)");
//            energy += ee;
//        }  
//            
////        System.out.println(")");
//        return energy;
//    }
    
    
    //evaluate energy of substructure from i -> pair[i]
    private static double componentEnergy(RNAPrimarySequence sequence, int[] pairList, int i) {
//        System.out.println(" Call component energy for the base-pair (" + i + ", " + pairList[i] + ")");

        /* calculate energy of substructure enclosed by (i,j) */
        double ee = 0, energy = 0;

        boolean isHairpin_Multiloop = false;

        int j, p, q;

        j = pairList[i];

        p = i;
        q = j;
        while (p < q) {
            while (pairList[++p] == -1);
            while (pairList[--q] == -1);

            /* MULTILOOP_A hairpin or multi-loop */
            if ((pairList[q] != p) || (p > q)) {
                isHairpin_Multiloop = true;
                break;
            }

            //a stack/bulge/interior loop found            
            ee = loopEnergy(sequence, i, j, p, q);

            energy += ee;

            i = p;
            j = q;
        }

        if (isHairpin_Multiloop) {
            /* p,q don't pair must have found hairpin or multiloop */
            if (p > q) {
                /* MULTILOOP_A hair pin */
                ee = hairpinEnergy(sequence, i, j);

                energy += ee;

                return energy;
            }

            /* (i,j) is exterior pair of multiloop */
//            System.out.print("  A multi-loop enclosed by (" + i + ", " + j + ")");
            int branch = 0; //number of branching helices
            int numUnpairs = p - i - 1; //number of unpairs
            while (p < j) {
                branch += 1;

                /* add up the contributions of the substructures of the multiloop */
                energy += componentEnergy(sequence, pairList, p);

                /* search for next base pair in multiloop */
                p = pairList[p];
                while (pairList[++p] == -1) {
                    numUnpairs++;
                }
            }
//            System.out.print(" has numUnpairs = " + numUnpairs + ", branch = " + branch);
            ee = multiloopEnergy(sequence, numUnpairs, branch);
//            System.out.println("  with multi-loop energy " + ee);

            energy += ee;
        }
        return energy;
    }

    public static double loopEnergy(RNAPrimarySequence sequence, int i, int j, int p, int q) {
        int n1 = p - i - 1;
        int n2 = j - q - 1;

        int nl, ns;

        double energy;

        if (n1 > n2) {
            nl = n1;
            ns = n2;
        } else {
            nl = n2;
            ns = n1;
        }

        if (nl == 0) {
            /* stack */
//            System.out.print("  A stack enclosed by (" + i + ", " + j + ")");
            energy = stackEnergy(sequence, i, j, p, q);
        } else if (ns == 0) {
            /* bulge */
//            System.out.print("  A bulge enclosed by (" + i + ", " + j + ")");
            energy = bulgeEnergy(sequence, i, j, p, q, nl);
        } else {
            /* internal loop */
//            System.out.print("  An internal loop enclosed by (" + i + ", " + j + ")");
            energy = internalLoopEnergy(sequence, i, j, p, q, nl, ns);
        }
//        System.out.println(" has energy " + energy + " (kcal)");
        return energy;
    }

    public static double stackEnergy(RNAPrimarySequence sequence, int i, int j, int p, int q) {
        char iNucleotide = sequence.getNucleotide(i);
        char jNucleotide = sequence.getNucleotide(j);

        char pNucleotide = sequence.getNucleotide(p);
        char qNucleotide = sequence.getNucleotide(q);

        return lookupStackEnergy(iNucleotide, jNucleotide, pNucleotide, qNucleotide);
    }

    public static double bulgeEnergy(RNAPrimarySequence sequence, int i, int j, int p, int q, int size) {
        double energy = 0.0;

        double ee;
        ee = (size <= ENERGY_CONSTANTS.MAXLOOP) ? (ENERGY_CONSTANTS.bulge[size]) : (ENERGY_CONSTANTS.bulge[ENERGY_CONSTANTS.MAXLOOP] + ENERGY_CONSTANTS.lxc37 * Math.log(size / (double) ENERGY_CONSTANTS.MAXLOOP));
//        System.out.print(" (initiation energy for bulge size " + size + " = " + ee + " (kcal)");
        energy += ee;

        if (size == 1) {
            ee = stackEnergy(sequence, i, j, p, q);
//            System.out.print(", stacking bonus = " + ee + " (kcal) ");
            energy += ee;
        }
        //AU/GU end penalty
        else {
            char iNucleotide = sequence.getNucleotide(i);
            char jNucleotide = sequence.getNucleotide(j);
            if (checkAUEnd(iNucleotide, jNucleotide)) {
                ee = ENERGY_CONSTANTS.AU_end;
//                System.out.print(", AU end penalty by ( i = " + i + " , j = " + j + " ) = " + ee + " (kcal) ");
                energy += ee;
            } else if (checkGUEnd(iNucleotide, jNucleotide)) {
                ee = ENERGY_CONSTANTS.GU_end;
//                System.out.print(", GU end penalty by ( i = " + i + " , j = " + j + " ) = " + ee + " (kcal) ");
                energy += ee;
            }

            char pNucleotide = sequence.getNucleotide(p);
            char qNucleotide = sequence.getNucleotide(q);
            if (checkAUEnd(pNucleotide, qNucleotide)) {
                ee = ENERGY_CONSTANTS.AU_end;
//                System.out.print(", AU end penalty by ( p = " + p + " , q = " + q + " ) = " + ee + " (kcal) ");
                energy += ee;
            } else if (checkGUEnd(pNucleotide, qNucleotide)) {
                ee = ENERGY_CONSTANTS.GU_end;
//                System.out.print(", GU end penalty by ( p = " + p + " , q = " + q + " ) = " + ee + " (kcal) ");
                energy += ee;
            }
        }
//        System.out.print(")");
        return energy;
    }

    public static double internalLoopEnergy(RNAPrimarySequence sequence, int i, int j, int p, int q, int big, int small) {
        double energy = 0.0;

        int size = big + small;
        int difference = big - small;

        char iNucleotide = sequence.getNucleotide(i);
        char jNucleotide = sequence.getNucleotide(j);

        char pNucleotide = sequence.getNucleotide(p);
        char qNucleotide = sequence.getNucleotide(q);

        if (size == 2) {
            char xNucleotide = sequence.getNucleotide(i + 1);
            char yNucleotide = sequence.getNucleotide(i + 1);
            energy = lookupInternalLoop11Energy(iNucleotide, jNucleotide, pNucleotide, qNucleotide, xNucleotide, yNucleotide);
        } else {
            double ee;

            ee = (size <= ENERGY_CONSTANTS.MAXLOOP) ? (ENERGY_CONSTANTS.internal_loop[size]) : (ENERGY_CONSTANTS.internal_loop[ENERGY_CONSTANTS.MAXLOOP] + ENERGY_CONSTANTS.lxc37 * Math.log((size) / (double) ENERGY_CONSTANTS.MAXLOOP));
//            System.out.print(" (initiation energy for internal loop size " + size + " = " + ee + " (kcal)");
            energy += ee;

            ee = difference * ENERGY_CONSTANTS.ASSYMETRY;
//            System.out.print(", assymetry term = " + ee + " (kcal)");
            energy += ee;

            if ((iNucleotide == NUCLEOTIDE_CODE.A && jNucleotide == NUCLEOTIDE_CODE.U)
                    || (iNucleotide == NUCLEOTIDE_CODE.G && jNucleotide == NUCLEOTIDE_CODE.U)
                    || (pNucleotide == NUCLEOTIDE_CODE.A && qNucleotide == NUCLEOTIDE_CODE.U)
                    || (pNucleotide == NUCLEOTIDE_CODE.G && qNucleotide == NUCLEOTIDE_CODE.U)) {
                ee = ENERGY_CONSTANTS.AU_GU_Closure_Internal;
//                System.out.print(", AU/GU closure = " + ee + " (kcal)");
                energy += ee;
            }

            if (small >= 2) {
                /* mismatch 1 */
                //AG or GA first mismatch
                if (   (sequence.getNucleotide(i + 1) == NUCLEOTIDE_CODE.A && sequence.getNucleotide(j - 1) == NUCLEOTIDE_CODE.G)
                    || (sequence.getNucleotide(i + 1) == NUCLEOTIDE_CODE.G && sequence.getNucleotide(j - 1) == NUCLEOTIDE_CODE.A)) {
                    ee = ENERGY_CONSTANTS.GA_AG_First_Mismatch_Internal;
//                    System.out.print(", AG/GA first mismatch bonus = " + ee + " (kcal)");
                    energy += ee;
                } //UU first mismatch
                else if (sequence.getNucleotide(i + 1) == NUCLEOTIDE_CODE.U && sequence.getNucleotide(j - 1) == NUCLEOTIDE_CODE.U) {
                    ee = ENERGY_CONSTANTS.UU_First_Mismatch_Internal;
//                    System.out.print(", AG/GA first mismatch bonus = " + ee + " (kcal)");
                    energy += ee;
                }

                /* mismatch 2 */
                //AG or GA first mismatch
                if (   (sequence.getNucleotide(p + 1) == NUCLEOTIDE_CODE.A && sequence.getNucleotide(q - 1) == NUCLEOTIDE_CODE.G)
                    || (sequence.getNucleotide(p + 1) == NUCLEOTIDE_CODE.G && sequence.getNucleotide(q - 1) == NUCLEOTIDE_CODE.A)) {
                    ee = ENERGY_CONSTANTS.GA_AG_First_Mismatch_Internal;
//                    System.out.print(", AG/GA first mismatch bonus = " + ee + " (kcal)");
                    energy += ee;
                } //UU first mismatch
                else if (sequence.getNucleotide(p + 1) == NUCLEOTIDE_CODE.U && sequence.getNucleotide(q - 1) == NUCLEOTIDE_CODE.U) {
                    ee = ENERGY_CONSTANTS.UU_First_Mismatch_Internal;
//                    System.out.print(", AG/GA first mismatch bonus = " + ee + " (kcal)");
                    energy += ee;
                }
            }
//            System.out.print(" )");

        }
        return energy;
    }

    public static double hairpinEnergy(RNAPrimarySequence sequence, int i, int j) {
//        System.out.print("  A hairpin enclosed by (" + i + ", " + j + ")");
        double energy = 0.0;

        double ee;
        //initialization
        int size = j - i - 1;

        ee = (size <= ENERGY_CONSTANTS.MAXLOOP) ? (ENERGY_CONSTANTS.hairpin[size]) : (ENERGY_CONSTANTS.hairpin[ENERGY_CONSTANTS.MAXLOOP] + ENERGY_CONSTANTS.lxc37 * Math.log(size / (double) ENERGY_CONSTANTS.MAXLOOP));
//        System.out.print(" (initiation energy for hairpin size " + size + " = " + ee + " (kcal)");
        energy += ee;

        //all c panalty
        boolean all_C_Loop = true;
        for (int runIndex = i + 1; runIndex <= j - 1; runIndex++) {
            if (sequence.getNucleotide(runIndex) != 'C') {
                all_C_Loop = false;
                break;
            }
        }

        if (all_C_Loop) {
//            System.out.print(", all C loop penalty = ");
            if (size == 3) {
                ee = ENERGY_CONSTANTS.ALL_C_3;
            } else {
                ee = (size * ENERGY_CONSTANTS.ALL_C_LOOP_A + ENERGY_CONSTANTS.ALL_C_LOOP_B);
            }
//            System.out.print(ee + " (kcal)");
            energy += ee;
        }

        //AU/GU end penalty
        char iNucleotide = sequence.getNucleotide(i);
        char jNucleotide = sequence.getNucleotide(j);

        //AU/GU end penalty
        if (checkAUEnd(iNucleotide, jNucleotide)) {
            ee = ENERGY_CONSTANTS.AU_end;
//            System.out.print(", AU end penalty = " + ee + " (kcal)");
            energy += ee;
        } else if (checkGUEnd(iNucleotide, jNucleotide)) {
            ee = ENERGY_CONSTANTS.GU_end;
//            System.out.print(", GU end penalty = " + ee + " (kcal)");
            energy += ee;
        }

        //bonus
        if (size > 3) {
            //GU closure
            if (iNucleotide == NUCLEOTIDE_CODE.G && jNucleotide == NUCLEOTIDE_CODE.U 
            && (i >= 2 && sequence.getNucleotide(i-1) == NUCLEOTIDE_CODE.G && sequence.getNucleotide(i-2) == NUCLEOTIDE_CODE.G) )
            {
                ee = ENERGY_CONSTANTS.GU_Closure_Hairpin;
//                System.out.print(", GU closure = " + ee + " (kcal)");
                energy += ee;                
            }

            char xNucleotide = sequence.getNucleotide(i + 1);
            char yNucleotide = sequence.getNucleotide(j - 1);

            //terminal mismatch
            ee = lookupTerminalMismatchEnergy(iNucleotide, jNucleotide, xNucleotide, yNucleotide);
//            System.out.print(", terminal mismatch = " + ee + " (kcal)");
            energy += ee;

            //UU or GA first mismatch
            if ((xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.U)
                    || (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.A)) {
                ee = ENERGY_CONSTANTS.UU_GA_First_Mismatch_Hairpin;
//                System.out.print(", UU/GA first mismatch bonus = " + ee + " (kcal)");
                energy += ee;
            }

            //tetraloop
            if (size == 4) {
                String xyzw = sequence.getSubsequence(i + 1, j - 1);

                ee = lookupTetraloopEnergy(iNucleotide, jNucleotide, xyzw);
//                System.out.print(", tetraloop = " + ee + " (kcal)");
                energy += ee;
            }
        }
//        System.out.println(" ) has energy " + energy + " (kcal)");
        return energy;
    }

    public static double multiloopEnergy(RNAPrimarySequence sequence, int numUnpairs, int branch) {
        double energy = ENERGY_CONSTANTS.MULTILOOP_A + (ENERGY_CONSTANTS.MULTILOOP_B * numUnpairs) + (ENERGY_CONSTANTS.MULTILOOP_C * branch);
        return energy;
    }

    private static boolean checkAUEnd(char iNucleotide, char jNucleotide) {
        if ((iNucleotide == NUCLEOTIDE_CODE.A && jNucleotide == NUCLEOTIDE_CODE.U)
                || (iNucleotide == NUCLEOTIDE_CODE.U && jNucleotide == NUCLEOTIDE_CODE.A)) {
            return true;
        }
        return false;
    }

    private static boolean checkGUEnd(char iNucleotide, char jNucleotide) {
        if ((iNucleotide == NUCLEOTIDE_CODE.G && jNucleotide == NUCLEOTIDE_CODE.U)
                || (iNucleotide == NUCLEOTIDE_CODE.U && jNucleotide == NUCLEOTIDE_CODE.G)) {
            return true;
        }
        return false;
    }

    private static double lookupDanglingEnergy(char iNucleotide, char jNucleotide, char xNucleotide, boolean isFivePrime) {
        //AU
        if (iNucleotide == NUCLEOTIDE_CODE.A && jNucleotide == NUCLEOTIDE_CODE.U) {
            if (xNucleotide == NUCLEOTIDE_CODE.A) {
                if (isFivePrime) {
                    return ENERGY_CONSTANTS.dangling_5[0][0];
                } else {
                    return ENERGY_CONSTANTS.dangling_3[0][0];
                }
            }
            else if (xNucleotide == NUCLEOTIDE_CODE.C) {
                if (isFivePrime) {
                    return ENERGY_CONSTANTS.dangling_5[0][1];
                } else {
                    return ENERGY_CONSTANTS.dangling_3[0][1];
                }
            }
            else if (xNucleotide == NUCLEOTIDE_CODE.G) {
                if (isFivePrime) {
                    return ENERGY_CONSTANTS.dangling_5[0][2];
                } else {
                    return ENERGY_CONSTANTS.dangling_3[0][2];
                }
            }
            else if (xNucleotide == NUCLEOTIDE_CODE.U) {
                if (isFivePrime) {
                    return ENERGY_CONSTANTS.dangling_5[0][3];
                } else {
                    return ENERGY_CONSTANTS.dangling_3[0][3];
                }
            }
        } //CG
        else if (iNucleotide == NUCLEOTIDE_CODE.C && jNucleotide == NUCLEOTIDE_CODE.G) {
            if (xNucleotide == NUCLEOTIDE_CODE.A) {
                if (isFivePrime) {
                    return ENERGY_CONSTANTS.dangling_5[1][0];
                } else {
                    return ENERGY_CONSTANTS.dangling_3[1][0];
                }
            }
            else if (xNucleotide == NUCLEOTIDE_CODE.C) {
                if (isFivePrime) {
                    return ENERGY_CONSTANTS.dangling_5[1][1];
                } else {
                    return ENERGY_CONSTANTS.dangling_3[1][1];
                }
            }
            else if (xNucleotide == NUCLEOTIDE_CODE.G) {
                if (isFivePrime) {
                    return ENERGY_CONSTANTS.dangling_5[1][2];
                } else {
                    return ENERGY_CONSTANTS.dangling_3[1][2];
                }
            }
            else if (xNucleotide == NUCLEOTIDE_CODE.U) {
                if (isFivePrime) {
                    return ENERGY_CONSTANTS.dangling_5[1][3];
                } else {
                    return ENERGY_CONSTANTS.dangling_3[1][3];
                }
            }
        } //GC
        else if (iNucleotide == NUCLEOTIDE_CODE.G && jNucleotide == NUCLEOTIDE_CODE.C) {
            if (xNucleotide == NUCLEOTIDE_CODE.A) {
                if (isFivePrime) {
                    return ENERGY_CONSTANTS.dangling_5[2][0];
                } else {
                    return ENERGY_CONSTANTS.dangling_3[2][0];
                }
            }
            else if (xNucleotide == NUCLEOTIDE_CODE.C) {
                if (isFivePrime) {
                    return ENERGY_CONSTANTS.dangling_5[2][1];
                } else {
                    return ENERGY_CONSTANTS.dangling_3[2][1];
                }
            }
            else if (xNucleotide == NUCLEOTIDE_CODE.G) {
                if (isFivePrime) {
                    return ENERGY_CONSTANTS.dangling_5[2][2];
                } else {
                    return ENERGY_CONSTANTS.dangling_3[2][2];
                }
            }
            else if (xNucleotide == NUCLEOTIDE_CODE.U) {
                if (isFivePrime) {
                    return ENERGY_CONSTANTS.dangling_5[2][3];
                } else {
                    return ENERGY_CONSTANTS.dangling_3[2][3];
                }
            }
        } //GU
        else if (iNucleotide == NUCLEOTIDE_CODE.G && jNucleotide == NUCLEOTIDE_CODE.U) {
            if (xNucleotide == NUCLEOTIDE_CODE.A) {
                if (isFivePrime) {
                    return ENERGY_CONSTANTS.dangling_5[3][0];
                } else {
                    return ENERGY_CONSTANTS.dangling_3[3][0];
                }
            }
            else if (xNucleotide == NUCLEOTIDE_CODE.C) {
                if (isFivePrime) {
                    return ENERGY_CONSTANTS.dangling_5[3][1];
                } else {
                    return ENERGY_CONSTANTS.dangling_3[3][1];
                }
            }
            else if (xNucleotide == NUCLEOTIDE_CODE.G) {
                if (isFivePrime) {
                    return ENERGY_CONSTANTS.dangling_5[3][2];
                } else {
                    return ENERGY_CONSTANTS.dangling_3[3][2];
                }
            }
            else if (xNucleotide == NUCLEOTIDE_CODE.U) {
                if (isFivePrime) {
                    return ENERGY_CONSTANTS.dangling_5[3][3];
                } else {
                    return ENERGY_CONSTANTS.dangling_3[3][3];
                }
            }
        } //UA
        else if (iNucleotide == NUCLEOTIDE_CODE.U && jNucleotide == NUCLEOTIDE_CODE.A) {
            if (xNucleotide == NUCLEOTIDE_CODE.A) {
                if (isFivePrime) {
                    return ENERGY_CONSTANTS.dangling_5[4][0];
                } else {
                    return ENERGY_CONSTANTS.dangling_3[4][0];
                }
            }
            else if (xNucleotide == NUCLEOTIDE_CODE.C) {
                if (isFivePrime) {
                    return ENERGY_CONSTANTS.dangling_5[4][1];
                } else {
                    return ENERGY_CONSTANTS.dangling_3[4][1];
                }
            }
            else if (xNucleotide == NUCLEOTIDE_CODE.G) {
                if (isFivePrime) {
                    return ENERGY_CONSTANTS.dangling_5[4][2];
                } else {
                    return ENERGY_CONSTANTS.dangling_3[4][2];
                }
            }
            else if (xNucleotide == NUCLEOTIDE_CODE.U) {
                if (isFivePrime) {
                    return ENERGY_CONSTANTS.dangling_5[4][3];
                } else {
                    return ENERGY_CONSTANTS.dangling_3[4][3];
                }
            }
        }//UG
        else if (iNucleotide == NUCLEOTIDE_CODE.U && jNucleotide == NUCLEOTIDE_CODE.G) {
            if (xNucleotide == NUCLEOTIDE_CODE.A) {
                if (isFivePrime) {
                    return ENERGY_CONSTANTS.dangling_5[5][0];
                } else {
                    return ENERGY_CONSTANTS.dangling_3[5][0];
                }
            }
            else if (xNucleotide == NUCLEOTIDE_CODE.C) {
                if (isFivePrime) {
                    return ENERGY_CONSTANTS.dangling_5[5][1];
                } else {
                    return ENERGY_CONSTANTS.dangling_3[5][1];
                }
            }
            else if (xNucleotide == NUCLEOTIDE_CODE.G) {
                if (isFivePrime) {
                    return ENERGY_CONSTANTS.dangling_5[5][2];
                } else {
                    return ENERGY_CONSTANTS.dangling_3[5][2];
                }
            }
            else if (xNucleotide == NUCLEOTIDE_CODE.U) {
                if (isFivePrime) {
                    return ENERGY_CONSTANTS.dangling_5[5][3];
                } else {
                    return ENERGY_CONSTANTS.dangling_3[5][3];
                }
            }
        }
        return 0.0;
    }

    private static double lookupTerminalMismatchEnergy(char iNucleotide, char jNucleotide, char xNucleotide, char yNucleotide) {
        //AU
        if (iNucleotide == NUCLEOTIDE_CODE.A && jNucleotide == NUCLEOTIDE_CODE.U) {
            //AA
            if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.terminal_mismatch[0][0];
            } //AC
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.terminal_mismatch[0][1];
            } //AG
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.terminal_mismatch[0][2];
            } //AU
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.terminal_mismatch[0][3];
            } //CA
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.terminal_mismatch[0][4];
            } //CC
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.terminal_mismatch[0][5];
            } //CG
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.terminal_mismatch[0][6];
            } //CU
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.terminal_mismatch[0][7];
            } //GA
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.terminal_mismatch[0][8];
            } //GC
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.terminal_mismatch[0][9];
            } //GG
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.terminal_mismatch[0][10];
            } //GU
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.terminal_mismatch[0][11];
            } //UA
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.terminal_mismatch[0][12];
            } //UC
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.terminal_mismatch[0][13];
            } //UG
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.terminal_mismatch[0][14];
            } //UU
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.terminal_mismatch[0][15];
            } else {
                throw new RuntimeException("Invalid code for terminal mismatch energy");
            }
        } //CG
        else if (iNucleotide == NUCLEOTIDE_CODE.C && jNucleotide == NUCLEOTIDE_CODE.G) {
            //AA
            if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.terminal_mismatch[1][0];
            } //AC
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.terminal_mismatch[1][1];
            } //AG
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.terminal_mismatch[1][2];
            } //AU
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.terminal_mismatch[1][3];
            } //CA
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.terminal_mismatch[1][4];
            } //CC
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.terminal_mismatch[1][5];
            } //CG
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.terminal_mismatch[1][6];
            } //CU
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.terminal_mismatch[1][7];
            } //GA
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.terminal_mismatch[1][8];
            } //GC
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.terminal_mismatch[1][9];
            } //GG
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.terminal_mismatch[1][10];
            } //GU
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.terminal_mismatch[1][11];
            } //UA
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.terminal_mismatch[1][12];
            } //UC
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.terminal_mismatch[1][13];
            } //UG
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.terminal_mismatch[1][14];
            } //UU
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.terminal_mismatch[1][15];
            } else {
                throw new RuntimeException("Invalid code for terminal mismatch energy");
            }
        } //GC
        else if (iNucleotide == NUCLEOTIDE_CODE.G && jNucleotide == NUCLEOTIDE_CODE.C) {
            //AA
            if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.terminal_mismatch[2][0];
            } //AC
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.terminal_mismatch[2][1];
            } //AG
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.terminal_mismatch[2][2];
            } //AU
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.terminal_mismatch[2][3];
            } //CA
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.terminal_mismatch[2][4];
            } //CC
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.terminal_mismatch[2][5];
            } //CG
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.terminal_mismatch[2][6];
            } //CU
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.terminal_mismatch[2][7];
            } //GA
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.terminal_mismatch[2][8];
            } //GC
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.terminal_mismatch[2][9];
            } //GG
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.terminal_mismatch[2][10];
            } //GU
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.terminal_mismatch[2][11];
            } //UA
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.terminal_mismatch[2][12];
            } //UC
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.terminal_mismatch[2][13];
            } //UG
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.terminal_mismatch[2][14];
            } //UU
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.terminal_mismatch[2][15];
            } else {
                throw new RuntimeException("Invalid code for terminal mismatch energy");
            }
        } //GU
        else if (iNucleotide == NUCLEOTIDE_CODE.G && jNucleotide == NUCLEOTIDE_CODE.U) {
            //AA
            if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.terminal_mismatch[3][0];
            } //AC
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.terminal_mismatch[3][1];
            } //AG
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.terminal_mismatch[3][2];
            } //AU
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.terminal_mismatch[3][3];
            } //CA
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.terminal_mismatch[3][4];
            } //CC
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.terminal_mismatch[3][5];
            } //CG
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.terminal_mismatch[3][6];
            } //CU
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.terminal_mismatch[3][7];
            } //GA
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.terminal_mismatch[3][8];
            } //GC
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.terminal_mismatch[3][9];
            } //GG
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.terminal_mismatch[3][10];
            } //GU
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.terminal_mismatch[3][11];
            } //UA
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.terminal_mismatch[3][12];
            } //UC
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.terminal_mismatch[3][13];
            } //UG
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.terminal_mismatch[3][14];
            } //UU
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.terminal_mismatch[3][15];
            } else {
                throw new RuntimeException("Invalid code for terminal mismatch energy");
            }
        } //UA
        else if (iNucleotide == NUCLEOTIDE_CODE.U && jNucleotide == NUCLEOTIDE_CODE.A) {
            //AA
            if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.terminal_mismatch[4][0];
            } //AC
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.terminal_mismatch[4][1];
            } //AG
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.terminal_mismatch[4][2];
            } //AU
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.terminal_mismatch[4][3];
            } //CA
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.terminal_mismatch[4][4];
            } //CC
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.terminal_mismatch[4][5];
            } //CG
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.terminal_mismatch[4][6];
            } //CU
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.terminal_mismatch[4][7];
            } //GA
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.terminal_mismatch[4][8];
            } //GC
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.terminal_mismatch[4][9];
            } //GG
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.terminal_mismatch[4][10];
            } //GU
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.terminal_mismatch[4][11];
            } //UA
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.terminal_mismatch[4][12];
            } //UC
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.terminal_mismatch[4][13];
            } //UG
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.terminal_mismatch[4][14];
            } //UU
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.terminal_mismatch[4][15];
            } else {
                throw new RuntimeException("Invalid code for terminal mismatch energy");
            }
        } //UG
        else if (iNucleotide == NUCLEOTIDE_CODE.U && jNucleotide == NUCLEOTIDE_CODE.G) {
            //AA
            if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.terminal_mismatch[5][0];
            } //AC
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.terminal_mismatch[5][1];
            } //AG
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.terminal_mismatch[5][2];
            } //AU
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.terminal_mismatch[5][3];
            } //CA
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.terminal_mismatch[5][4];
            } //CC
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.terminal_mismatch[5][5];
            } //CG
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.terminal_mismatch[5][6];
            } //CU
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.terminal_mismatch[5][7];
            } //GA
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.terminal_mismatch[5][8];
            } //GC
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.terminal_mismatch[5][9];
            } //GG
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.terminal_mismatch[5][10];
            } //GU
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.terminal_mismatch[5][11];
            } //UA
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.terminal_mismatch[5][12];
            } //UC
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.terminal_mismatch[5][13];
            } //UG
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.terminal_mismatch[5][14];
            } //UU
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.terminal_mismatch[5][15];
            } else {
                throw new RuntimeException("Invalid code for terminal mismatch energy");
            }
        } else {
            throw new RuntimeException("Invalid code for terminal mismatch energy");
        }
    }

    private static double lookupTetraloopEnergy(char iNucleotide, char jNucleotide, String xyzw) {
        //GGGGAC
        if (iNucleotide == NUCLEOTIDE_CODE.G && jNucleotide == NUCLEOTIDE_CODE.C && xyzw.equals("GGGA")) {
            return ENERGY_CONSTANTS.tetraloop_bonus[0];
        } //GGUGAC
        else if (iNucleotide == NUCLEOTIDE_CODE.G && jNucleotide == NUCLEOTIDE_CODE.C && xyzw.equals("GUGA")) {
            return ENERGY_CONSTANTS.tetraloop_bonus[1];
        } //CGAAAG
        else if (iNucleotide == NUCLEOTIDE_CODE.C && jNucleotide == NUCLEOTIDE_CODE.G && xyzw.equals("GAAA")) {
            return ENERGY_CONSTANTS.tetraloop_bonus[2];
        } //GGAGAC
        else if (iNucleotide == NUCLEOTIDE_CODE.G && jNucleotide == NUCLEOTIDE_CODE.C && xyzw.equals("GAGA")) {
            return ENERGY_CONSTANTS.tetraloop_bonus[3];
        } //CGCAAG
        else if (iNucleotide == NUCLEOTIDE_CODE.C && jNucleotide == NUCLEOTIDE_CODE.G && xyzw.equals("GCAA")) {
            return ENERGY_CONSTANTS.tetraloop_bonus[4];
        } //GGAAAC
        else if (iNucleotide == NUCLEOTIDE_CODE.G && jNucleotide == NUCLEOTIDE_CODE.C && xyzw.equals("GAAA")) {
            return ENERGY_CONSTANTS.tetraloop_bonus[5];
        } //CGGAAG	
        else if (iNucleotide == NUCLEOTIDE_CODE.C && jNucleotide == NUCLEOTIDE_CODE.G && xyzw.equals("GGAA")) {
            return ENERGY_CONSTANTS.tetraloop_bonus[6];
        } //CUUCGG
        else if (iNucleotide == NUCLEOTIDE_CODE.C && jNucleotide == NUCLEOTIDE_CODE.G && xyzw.equals("UUCG")) {
            return ENERGY_CONSTANTS.tetraloop_bonus[7];
        } //CGUGAG
        else if (iNucleotide == NUCLEOTIDE_CODE.C && jNucleotide == NUCLEOTIDE_CODE.G && xyzw.equals("GUGA")) {
            return ENERGY_CONSTANTS.tetraloop_bonus[8];
        } //CGAAGG
        else if (iNucleotide == NUCLEOTIDE_CODE.C && jNucleotide == NUCLEOTIDE_CODE.G && xyzw.equals("GAAG")) {
            return ENERGY_CONSTANTS.tetraloop_bonus[9];
        } //CUACGG
        else if (iNucleotide == NUCLEOTIDE_CODE.C && jNucleotide == NUCLEOTIDE_CODE.G && xyzw.equals("UACG")) {
            return ENERGY_CONSTANTS.tetraloop_bonus[10];
        } //GGCAAC
        else if (iNucleotide == NUCLEOTIDE_CODE.G && jNucleotide == NUCLEOTIDE_CODE.C && xyzw.equals("GCAA")) {
            return ENERGY_CONSTANTS.tetraloop_bonus[11];
        } //CGCGAG
        else if (iNucleotide == NUCLEOTIDE_CODE.C && jNucleotide == NUCLEOTIDE_CODE.G && xyzw.equals("GCGA")) {
            return ENERGY_CONSTANTS.tetraloop_bonus[12];
        } //UGAGAG
        else if (iNucleotide == NUCLEOTIDE_CODE.U && jNucleotide == NUCLEOTIDE_CODE.G && xyzw.equals("GAGA")) {
            return ENERGY_CONSTANTS.tetraloop_bonus[13];
        } //CGAGAG
        else if (iNucleotide == NUCLEOTIDE_CODE.C && jNucleotide == NUCLEOTIDE_CODE.G && xyzw.equals("GAGA")) {
            return ENERGY_CONSTANTS.tetraloop_bonus[14];
        } //AGAAAU
        else if (iNucleotide == NUCLEOTIDE_CODE.A && jNucleotide == NUCLEOTIDE_CODE.U && xyzw.equals("GAAA")) {
            return ENERGY_CONSTANTS.tetraloop_bonus[15];
        } //CGUAAG
        else if (iNucleotide == NUCLEOTIDE_CODE.C && jNucleotide == NUCLEOTIDE_CODE.G && xyzw.equals("GUAA")) {
            return ENERGY_CONSTANTS.tetraloop_bonus[16];
        } //CUAACG
        else if (iNucleotide == NUCLEOTIDE_CODE.C && jNucleotide == NUCLEOTIDE_CODE.G && xyzw.equals("UAAC")) {
            return ENERGY_CONSTANTS.tetraloop_bonus[17];
        } //UGAAAG
        else if (iNucleotide == NUCLEOTIDE_CODE.U && jNucleotide == NUCLEOTIDE_CODE.G && xyzw.equals("GAAA")) {
            return ENERGY_CONSTANTS.tetraloop_bonus[18];
        } //GGAAGC
        else if (iNucleotide == NUCLEOTIDE_CODE.G && jNucleotide == NUCLEOTIDE_CODE.C && xyzw.equals("GAAG")) {
            return ENERGY_CONSTANTS.tetraloop_bonus[19];
        } //GGGAAC
        else if (iNucleotide == NUCLEOTIDE_CODE.G && jNucleotide == NUCLEOTIDE_CODE.C && xyzw.equals("GGAA")) {
            return ENERGY_CONSTANTS.tetraloop_bonus[20];
        } //UGAAAA
        else if (iNucleotide == NUCLEOTIDE_CODE.U && jNucleotide == NUCLEOTIDE_CODE.A && xyzw.equals("GAAA")) {
            return ENERGY_CONSTANTS.tetraloop_bonus[21];
        } //AGCAAU
        else if (iNucleotide == NUCLEOTIDE_CODE.A && jNucleotide == NUCLEOTIDE_CODE.U && xyzw.equals("GCAA")) {
            return ENERGY_CONSTANTS.tetraloop_bonus[22];
        } //AGUAAU
        else if (iNucleotide == NUCLEOTIDE_CODE.A && jNucleotide == NUCLEOTIDE_CODE.U && xyzw.equals("GUAA")) {
            return ENERGY_CONSTANTS.tetraloop_bonus[23];
        } //CGGGAG
        else if (iNucleotide == NUCLEOTIDE_CODE.C && jNucleotide == NUCLEOTIDE_CODE.G && xyzw.equals("GGGA")) {
            return ENERGY_CONSTANTS.tetraloop_bonus[24];
        } //AGUGAU
        else if (iNucleotide == NUCLEOTIDE_CODE.A && jNucleotide == NUCLEOTIDE_CODE.U && xyzw.equals("GUGA")) {
            return ENERGY_CONSTANTS.tetraloop_bonus[25];
        } //GGCGAC
        else if (iNucleotide == NUCLEOTIDE_CODE.G && jNucleotide == NUCLEOTIDE_CODE.C && xyzw.equals("GCGA")) {
            return ENERGY_CONSTANTS.tetraloop_bonus[26];
        } //GGGAGC
        else if (iNucleotide == NUCLEOTIDE_CODE.G && jNucleotide == NUCLEOTIDE_CODE.C && xyzw.equals("GGAG")) {
            return ENERGY_CONSTANTS.tetraloop_bonus[27];
        } //GUGAAC
        else if (iNucleotide == NUCLEOTIDE_CODE.G && jNucleotide == NUCLEOTIDE_CODE.C && xyzw.equals("UGAA")) {
            return ENERGY_CONSTANTS.tetraloop_bonus[28];
        } //UGGAAA
        else if (iNucleotide == NUCLEOTIDE_CODE.U && jNucleotide == NUCLEOTIDE_CODE.A && xyzw.equals("GGAA")) {
            return ENERGY_CONSTANTS.tetraloop_bonus[29];
        } else {
            return 0;
        }
    }

    private static double lookupStackEnergy(char iNucleotide, char jNucleotide, char pNucleotide, char qNucleotide) {
        //CG
        if ((iNucleotide == NUCLEOTIDE_CODE.C && jNucleotide == NUCLEOTIDE_CODE.G)
                && (pNucleotide == NUCLEOTIDE_CODE.C && qNucleotide == NUCLEOTIDE_CODE.G)) {
            return ENERGY_CONSTANTS.stack[0][0];
        } else if ((iNucleotide == NUCLEOTIDE_CODE.C && jNucleotide == NUCLEOTIDE_CODE.G)
                && (pNucleotide == NUCLEOTIDE_CODE.G && qNucleotide == NUCLEOTIDE_CODE.C)) {
            return ENERGY_CONSTANTS.stack[0][1];
        } else if ((iNucleotide == NUCLEOTIDE_CODE.C && jNucleotide == NUCLEOTIDE_CODE.G)
                && (pNucleotide == NUCLEOTIDE_CODE.G && qNucleotide == NUCLEOTIDE_CODE.U)) {
            return ENERGY_CONSTANTS.stack[0][2];
        } else if ((iNucleotide == NUCLEOTIDE_CODE.C && jNucleotide == NUCLEOTIDE_CODE.G)
                && (pNucleotide == NUCLEOTIDE_CODE.U && qNucleotide == NUCLEOTIDE_CODE.G)) {
            return ENERGY_CONSTANTS.stack[0][3];
        } else if ((iNucleotide == NUCLEOTIDE_CODE.C && jNucleotide == NUCLEOTIDE_CODE.G)
                && (pNucleotide == NUCLEOTIDE_CODE.A && qNucleotide == NUCLEOTIDE_CODE.U)) {
            return ENERGY_CONSTANTS.stack[0][4];
        } else if ((iNucleotide == NUCLEOTIDE_CODE.C && jNucleotide == NUCLEOTIDE_CODE.G)
                && (pNucleotide == NUCLEOTIDE_CODE.U && qNucleotide == NUCLEOTIDE_CODE.A)) {
            return ENERGY_CONSTANTS.stack[0][5];
        } //GC
        else if ((iNucleotide == NUCLEOTIDE_CODE.G && jNucleotide == NUCLEOTIDE_CODE.C)
                && (pNucleotide == NUCLEOTIDE_CODE.C && qNucleotide == NUCLEOTIDE_CODE.G)) {
            return ENERGY_CONSTANTS.stack[1][0];
        } else if ((iNucleotide == NUCLEOTIDE_CODE.G && jNucleotide == NUCLEOTIDE_CODE.C)
                && (pNucleotide == NUCLEOTIDE_CODE.G && qNucleotide == NUCLEOTIDE_CODE.C)) {
            return ENERGY_CONSTANTS.stack[1][1];
        } else if ((iNucleotide == NUCLEOTIDE_CODE.G && jNucleotide == NUCLEOTIDE_CODE.C)
                && (pNucleotide == NUCLEOTIDE_CODE.G && qNucleotide == NUCLEOTIDE_CODE.U)) {
            return ENERGY_CONSTANTS.stack[1][2];
        } else if ((iNucleotide == NUCLEOTIDE_CODE.G && jNucleotide == NUCLEOTIDE_CODE.C)
                && (pNucleotide == NUCLEOTIDE_CODE.U && qNucleotide == NUCLEOTIDE_CODE.G)) {
            return ENERGY_CONSTANTS.stack[1][3];
        } else if ((iNucleotide == NUCLEOTIDE_CODE.G && jNucleotide == NUCLEOTIDE_CODE.C)
                && (pNucleotide == NUCLEOTIDE_CODE.A && qNucleotide == NUCLEOTIDE_CODE.U)) {
            return ENERGY_CONSTANTS.stack[1][4];
        } else if ((iNucleotide == NUCLEOTIDE_CODE.G && jNucleotide == NUCLEOTIDE_CODE.C)
                && (pNucleotide == NUCLEOTIDE_CODE.U && qNucleotide == NUCLEOTIDE_CODE.A)) {
            return ENERGY_CONSTANTS.stack[1][5];
        } //GU
        else if ((iNucleotide == NUCLEOTIDE_CODE.G && jNucleotide == NUCLEOTIDE_CODE.U)
                && (pNucleotide == NUCLEOTIDE_CODE.C && qNucleotide == NUCLEOTIDE_CODE.G)) {
            return ENERGY_CONSTANTS.stack[2][0];
        } else if ((iNucleotide == NUCLEOTIDE_CODE.G && jNucleotide == NUCLEOTIDE_CODE.U)
                && (pNucleotide == NUCLEOTIDE_CODE.G && qNucleotide == NUCLEOTIDE_CODE.C)) {
            return ENERGY_CONSTANTS.stack[2][1];
        } else if ((iNucleotide == NUCLEOTIDE_CODE.G && jNucleotide == NUCLEOTIDE_CODE.U)
                && (pNucleotide == NUCLEOTIDE_CODE.G && qNucleotide == NUCLEOTIDE_CODE.U)) {
            return ENERGY_CONSTANTS.stack[2][2];
        } else if ((iNucleotide == NUCLEOTIDE_CODE.G && jNucleotide == NUCLEOTIDE_CODE.U)
                && (pNucleotide == NUCLEOTIDE_CODE.U && qNucleotide == NUCLEOTIDE_CODE.G)) {
            return ENERGY_CONSTANTS.stack[2][3];
        } else if ((iNucleotide == NUCLEOTIDE_CODE.G && jNucleotide == NUCLEOTIDE_CODE.U)
                && (pNucleotide == NUCLEOTIDE_CODE.A && qNucleotide == NUCLEOTIDE_CODE.U)) {
            return ENERGY_CONSTANTS.stack[2][4];
        } else if ((iNucleotide == NUCLEOTIDE_CODE.G && jNucleotide == NUCLEOTIDE_CODE.U)
                && (pNucleotide == NUCLEOTIDE_CODE.U && qNucleotide == NUCLEOTIDE_CODE.A)) {
            return ENERGY_CONSTANTS.stack[2][5];
        } //UG
        else if ((iNucleotide == NUCLEOTIDE_CODE.U && jNucleotide == NUCLEOTIDE_CODE.G)
                && (pNucleotide == NUCLEOTIDE_CODE.C && qNucleotide == NUCLEOTIDE_CODE.G)) {
            return ENERGY_CONSTANTS.stack[3][0];
        } else if ((iNucleotide == NUCLEOTIDE_CODE.U && jNucleotide == NUCLEOTIDE_CODE.G)
                && (pNucleotide == NUCLEOTIDE_CODE.G && qNucleotide == NUCLEOTIDE_CODE.C)) {
            return ENERGY_CONSTANTS.stack[3][1];
        } else if ((iNucleotide == NUCLEOTIDE_CODE.U && jNucleotide == NUCLEOTIDE_CODE.G)
                && (pNucleotide == NUCLEOTIDE_CODE.G && qNucleotide == NUCLEOTIDE_CODE.U)) {
            return ENERGY_CONSTANTS.stack[3][2];
        } else if ((iNucleotide == NUCLEOTIDE_CODE.U && jNucleotide == NUCLEOTIDE_CODE.G)
                && (pNucleotide == NUCLEOTIDE_CODE.U && qNucleotide == NUCLEOTIDE_CODE.G)) {
            return ENERGY_CONSTANTS.stack[3][3];
        } else if ((iNucleotide == NUCLEOTIDE_CODE.U && jNucleotide == NUCLEOTIDE_CODE.G)
                && (pNucleotide == NUCLEOTIDE_CODE.A && qNucleotide == NUCLEOTIDE_CODE.U)) {
            return ENERGY_CONSTANTS.stack[3][4];
        } else if ((iNucleotide == NUCLEOTIDE_CODE.U && jNucleotide == NUCLEOTIDE_CODE.G)
                && (pNucleotide == NUCLEOTIDE_CODE.U && qNucleotide == NUCLEOTIDE_CODE.A)) {
            return ENERGY_CONSTANTS.stack[3][5];
        } //AU
        else if ((iNucleotide == NUCLEOTIDE_CODE.A && jNucleotide == NUCLEOTIDE_CODE.U)
                && (pNucleotide == NUCLEOTIDE_CODE.C && qNucleotide == NUCLEOTIDE_CODE.G)) {
            return ENERGY_CONSTANTS.stack[4][0];
        } else if ((iNucleotide == NUCLEOTIDE_CODE.A && jNucleotide == NUCLEOTIDE_CODE.U)
                && (pNucleotide == NUCLEOTIDE_CODE.G && qNucleotide == NUCLEOTIDE_CODE.C)) {
            return ENERGY_CONSTANTS.stack[4][1];
        } else if ((iNucleotide == NUCLEOTIDE_CODE.A && jNucleotide == NUCLEOTIDE_CODE.U)
                && (pNucleotide == NUCLEOTIDE_CODE.G && qNucleotide == NUCLEOTIDE_CODE.U)) {
            return ENERGY_CONSTANTS.stack[4][2];
        } else if ((iNucleotide == NUCLEOTIDE_CODE.A && jNucleotide == NUCLEOTIDE_CODE.U)
                && (pNucleotide == NUCLEOTIDE_CODE.U && qNucleotide == NUCLEOTIDE_CODE.G)) {
            return ENERGY_CONSTANTS.stack[4][3];
        } else if ((iNucleotide == NUCLEOTIDE_CODE.A && jNucleotide == NUCLEOTIDE_CODE.U)
                && (pNucleotide == NUCLEOTIDE_CODE.A && qNucleotide == NUCLEOTIDE_CODE.U)) {
            return ENERGY_CONSTANTS.stack[4][4];
        } else if ((iNucleotide == NUCLEOTIDE_CODE.A && jNucleotide == NUCLEOTIDE_CODE.U)
                && (pNucleotide == NUCLEOTIDE_CODE.U && qNucleotide == NUCLEOTIDE_CODE.A)) {
            return ENERGY_CONSTANTS.stack[4][5];
        } //UA
        else if ((iNucleotide == NUCLEOTIDE_CODE.U && jNucleotide == NUCLEOTIDE_CODE.A)
                && (pNucleotide == NUCLEOTIDE_CODE.C && qNucleotide == NUCLEOTIDE_CODE.G)) {
            return ENERGY_CONSTANTS.stack[5][0];
        } else if ((iNucleotide == NUCLEOTIDE_CODE.U && jNucleotide == NUCLEOTIDE_CODE.A)
                && (pNucleotide == NUCLEOTIDE_CODE.G && qNucleotide == NUCLEOTIDE_CODE.C)) {
            return ENERGY_CONSTANTS.stack[5][1];
        } else if ((iNucleotide == NUCLEOTIDE_CODE.U && jNucleotide == NUCLEOTIDE_CODE.A)
                && (pNucleotide == NUCLEOTIDE_CODE.G && qNucleotide == NUCLEOTIDE_CODE.U)) {
            return ENERGY_CONSTANTS.stack[5][2];
        } else if ((iNucleotide == NUCLEOTIDE_CODE.U && jNucleotide == NUCLEOTIDE_CODE.A)
                && (pNucleotide == NUCLEOTIDE_CODE.U && qNucleotide == NUCLEOTIDE_CODE.G)) {
            return ENERGY_CONSTANTS.stack[5][3];
        } else if ((iNucleotide == NUCLEOTIDE_CODE.U && jNucleotide == NUCLEOTIDE_CODE.A)
                && (pNucleotide == NUCLEOTIDE_CODE.A && qNucleotide == NUCLEOTIDE_CODE.U)) {
            return ENERGY_CONSTANTS.stack[5][4];
        } else if ((iNucleotide == NUCLEOTIDE_CODE.U && jNucleotide == NUCLEOTIDE_CODE.A)
                && (pNucleotide == NUCLEOTIDE_CODE.U && qNucleotide == NUCLEOTIDE_CODE.A)) {
            return ENERGY_CONSTANTS.stack[5][5];
        } else {
            throw new RuntimeException("Invalid code for stack energy");
        }
    }

    private static double lookupInternalLoop11Energy(char iNucleotide, char jNucleotide, char pNucleotide, char qNucleotide, char xNucleotide, char yNucleotide) {
        //AU-
        /*AU-AU*/
        if (iNucleotide == NUCLEOTIDE_CODE.A && jNucleotide == NUCLEOTIDE_CODE.U && pNucleotide == NUCLEOTIDE_CODE.A && qNucleotide == NUCLEOTIDE_CODE.U) {
            //AA
            if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[0][0];
            } //AC
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[0][1];
            } //AG
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[0][2];
            } //AU
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[0][3];
            } //CA
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[0][4];
            } //CC
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[0][5];
            } //CG
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[0][6];
            } //CU
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[0][7];
            } //GA
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[0][8];
            } //GC
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[0][9];
            } //GG
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[0][10];
            } //GU
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[0][11];
            } //UA
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[0][12];
            } //UC
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[0][13];
            } //UG
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[0][14];
            } //UU
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[0][15];
            } else {
                throw new RuntimeException("Invalid code for internal loop 1x1");
            }
        } /*AU-CG*/ else if (iNucleotide == NUCLEOTIDE_CODE.A && jNucleotide == NUCLEOTIDE_CODE.U && pNucleotide == NUCLEOTIDE_CODE.C && qNucleotide == NUCLEOTIDE_CODE.G) {
            //AA
            if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[1][0];
            } //AC
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[1][1];
            } //AG
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[1][2];
            } //AU
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[1][3];
            } //CA
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[1][4];
            } //CC
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[1][5];
            } //CG
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[1][6];
            } //CU
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[1][7];
            } //GA
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[1][8];
            } //GC
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[1][9];
            } //GG
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[1][10];
            } //GU
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[1][11];
            } //UA
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[1][12];
            } //UC
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[1][13];
            } //UG
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[1][14];
            } //UU
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[1][15];
            } else {
                throw new RuntimeException("Invalid code for internal loop 1x1");
            }
        } /*AU-GC*/ else if (iNucleotide == NUCLEOTIDE_CODE.A && jNucleotide == NUCLEOTIDE_CODE.U && pNucleotide == NUCLEOTIDE_CODE.G && qNucleotide == NUCLEOTIDE_CODE.C) {
            //AA
            if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[2][0];
            } //AC
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[2][1];
            } //AG
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[2][2];
            } //AU
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[2][3];
            } //CA
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[2][4];
            } //CC
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[2][5];
            } //CG
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[2][6];
            } //CU
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[2][7];
            } //GA
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[2][8];
            } //GC
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[2][9];
            } //GG
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[2][10];
            } //GU
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[2][11];
            } //UA
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[2][12];
            } //UC
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[2][13];
            } //UG
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[2][14];
            } //UU
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[2][15];
            } else {
                throw new RuntimeException("Invalid code for internal loop 1x1");
            }
        } /*AU-UA*/ else if (iNucleotide == NUCLEOTIDE_CODE.A && jNucleotide == NUCLEOTIDE_CODE.U && pNucleotide == NUCLEOTIDE_CODE.U && qNucleotide == NUCLEOTIDE_CODE.A) {
            //AA
            if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[3][0];
            } //AC
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[3][1];
            } //AG
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[3][2];
            } //AU
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[3][3];
            } //CA
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[3][4];
            } //CC
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[3][5];
            } //CG
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[3][6];
            } //CU
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[3][7];
            } //GA
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[3][8];
            } //GC
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[3][9];
            } //GG
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[3][10];
            } //GU
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[3][11];
            } //UA
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[3][12];
            } //UC
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[3][13];
            } //UG
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[3][14];
            } //UU
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[3][15];
            } else {
                throw new RuntimeException("Invalid code for internal loop 1x1");
            }
        } /*AU-GU*/ else if (iNucleotide == NUCLEOTIDE_CODE.A && jNucleotide == NUCLEOTIDE_CODE.U && pNucleotide == NUCLEOTIDE_CODE.G && qNucleotide == NUCLEOTIDE_CODE.U) {
            //AA
            if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[4][0];
            } //AC
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[4][1];
            } //AG
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[4][2];
            } //AU
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[4][3];
            } //CA
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[4][4];
            } //CC
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[4][5];
            } //CG
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[4][6];
            } //CU
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[4][7];
            } //GA
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[4][8];
            } //GC
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[4][9];
            } //GG
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[4][10];
            } //GU
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[4][11];
            } //UA
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[4][12];
            } //UC
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[4][13];
            } //UG
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[4][14];
            } //UU
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[4][15];
            } else {
                throw new RuntimeException("Invalid code for internal loop 1x1");
            }
        } /*AU-UG*/ else if (iNucleotide == NUCLEOTIDE_CODE.A && jNucleotide == NUCLEOTIDE_CODE.U && pNucleotide == NUCLEOTIDE_CODE.U && qNucleotide == NUCLEOTIDE_CODE.G) {
            //AA
            if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[5][0];
            } //AC
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[5][1];
            } //AG
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[5][2];
            } //AU
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[5][3];
            } //CA
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[5][4];
            } //CC
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[5][5];
            } //CG
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[5][6];
            } //CU
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[5][7];
            } //GA
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[5][8];
            } //GC
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[5][9];
            } //GG
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[5][10];
            } //GU
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[5][11];
            } //UA
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[5][12];
            } //UC
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[5][13];
            } //UG
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[5][14];
            } //UU
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[5][15];
            } else {
                throw new RuntimeException("Invalid code for internal loop 1x1");
            }
        } //CG-
        /*CG-AU*/ else if (iNucleotide == NUCLEOTIDE_CODE.C && jNucleotide == NUCLEOTIDE_CODE.G && pNucleotide == NUCLEOTIDE_CODE.A && qNucleotide == NUCLEOTIDE_CODE.U) {
            //AA
            if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[6][0];
            } //AC
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[6][1];
            } //AG
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[6][2];
            } //AU
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[6][3];
            } //CA
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[6][4];
            } //CC
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[6][5];
            } //CG
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[6][6];
            } //CU
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[6][7];
            } //GA
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[6][8];
            } //GC
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[6][9];
            } //GG
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[6][10];
            } //GU
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[6][11];
            } //UA
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[6][12];
            } //UC
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[6][13];
            } //UG
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[6][14];
            } //UU
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[6][15];
            } else {
                throw new RuntimeException("Invalid code for internal loop 1x1");
            }
        } /*CG-CG*/ else if (iNucleotide == NUCLEOTIDE_CODE.C && jNucleotide == NUCLEOTIDE_CODE.G && pNucleotide == NUCLEOTIDE_CODE.C && qNucleotide == NUCLEOTIDE_CODE.G) {
            //AA
            if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[7][0];
            } //AC
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[7][1];
            } //AG
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[7][2];
            } //AU
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[7][3];
            } //CA
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[7][4];
            } //CC
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[7][5];
            } //CG
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[7][6];
            } //CU
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[7][7];
            } //GA
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[7][8];
            } //GC
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[7][9];
            } //GG
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[7][10];
            } //GU
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[7][11];
            } //UA
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[7][12];
            } //UC
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[7][13];
            } //UG
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[7][14];
            } //UU
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[7][15];
            } else {
                throw new RuntimeException("Invalid code for internal loop 1x1");
            }
        } /*CG-GC*/ else if (iNucleotide == NUCLEOTIDE_CODE.C && jNucleotide == NUCLEOTIDE_CODE.G && pNucleotide == NUCLEOTIDE_CODE.G && qNucleotide == NUCLEOTIDE_CODE.C) {
            //AA
            if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[8][0];
            } //AC
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[8][1];
            } //AG
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[8][2];
            } //AU
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[8][3];
            } //CA
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[8][4];
            } //CC
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[8][5];
            } //CG
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[8][6];
            } //CU
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[8][7];
            } //GA
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[8][8];
            } //GC
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[8][9];
            } //GG
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[8][10];
            } //GU
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[8][11];
            } //UA
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[8][12];
            } //UC
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[8][13];
            } //UG
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[8][14];
            } //UU
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[8][15];
            } else {
                throw new RuntimeException("Invalid code for internal loop 1x1");
            }
        } /*CG-UA*/ else if (iNucleotide == NUCLEOTIDE_CODE.C && jNucleotide == NUCLEOTIDE_CODE.G && pNucleotide == NUCLEOTIDE_CODE.U && qNucleotide == NUCLEOTIDE_CODE.A) {
            //AA
            if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[9][0];
            } //AC
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[9][1];
            } //AG
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[9][2];
            } //AU
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[9][3];
            } //CA
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[9][4];
            } //CC
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[9][5];
            } //CG
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[9][6];
            } //CU
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[9][7];
            } //GA
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[9][8];
            } //GC
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[9][9];
            } //GG
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[9][10];
            } //GU
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[9][11];
            } //UA
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[9][12];
            } //UC
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[9][13];
            } //UG
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[9][14];
            } //UU
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[9][15];
            } else {
                throw new RuntimeException("Invalid code for internal loop 1x1");
            }
        } /*CG-GU*/ else if (iNucleotide == NUCLEOTIDE_CODE.C && jNucleotide == NUCLEOTIDE_CODE.G && pNucleotide == NUCLEOTIDE_CODE.G && qNucleotide == NUCLEOTIDE_CODE.U) {
            //AA
            if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[10][0];
            } //AC
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[10][1];
            } //AG
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[10][2];
            } //AU
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[10][3];
            } //CA
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[10][4];
            } //CC
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[10][5];
            } //CG
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[10][6];
            } //CU
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[10][7];
            } //GA
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[10][8];
            } //GC
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[10][9];
            } //GG
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[10][10];
            } //GU
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[10][11];
            } //UA
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[10][12];
            } //UC
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[10][13];
            } //UG
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[10][14];
            } //UU
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[10][15];
            } else {
                throw new RuntimeException("Invalid code for internal loop 1x1");
            }
        } /*CG-UG*/ else if (iNucleotide == NUCLEOTIDE_CODE.C && jNucleotide == NUCLEOTIDE_CODE.G && pNucleotide == NUCLEOTIDE_CODE.U && qNucleotide == NUCLEOTIDE_CODE.G) {
            //AA
            if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[11][0];
            } //AC
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[11][1];
            } //AG
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[11][2];
            } //AU
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[11][3];
            } //CA
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[11][4];
            } //CC
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[11][5];
            } //CG
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[11][6];
            } //CU
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[11][7];
            } //GA
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[11][8];
            } //GC
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[11][9];
            } //GG
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[11][10];
            } //GU
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[11][11];
            } //UA
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[11][12];
            } //UC
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[11][13];
            } //UG
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[11][14];
            } //UU
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[11][15];
            } else {
                throw new RuntimeException("Invalid code for internal loop 1x1");
            }
        } //GC
        /*GC-AU*/ else if (iNucleotide == NUCLEOTIDE_CODE.G && jNucleotide == NUCLEOTIDE_CODE.C && pNucleotide == NUCLEOTIDE_CODE.A && qNucleotide == NUCLEOTIDE_CODE.U) {
            //AA
            if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[12][0];
            } //AC
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[12][1];
            } //AG
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[12][2];
            } //AU
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[12][3];
            } //CA
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[12][4];
            } //CC
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[12][5];
            } //CG
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[12][6];
            } //CU
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[12][7];
            } //GA
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[12][8];
            } //GC
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[12][9];
            } //GG
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[12][10];
            } //GU
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[12][11];
            } //UA
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[12][12];
            } //UC
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[12][13];
            } //UG
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[12][14];
            } //UU
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[12][15];
            } else {
                throw new RuntimeException("Invalid code for internal loop 1x1");
            }
        } /*GC-CG*/ else if (iNucleotide == NUCLEOTIDE_CODE.G && jNucleotide == NUCLEOTIDE_CODE.C && pNucleotide == NUCLEOTIDE_CODE.C && qNucleotide == NUCLEOTIDE_CODE.G) {
            //AA
            if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[13][0];
            } //AC
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[13][1];
            } //AG
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[13][2];
            } //AU
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[13][3];
            } //CA
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[13][4];
            } //CC
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[13][5];
            } //CG
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[13][6];
            } //CU
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[13][7];
            } //GA
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[13][8];
            } //GC
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[13][9];
            } //GG
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[13][10];
            } //GU
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[13][11];
            } //UA
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[13][12];
            } //UC
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[13][13];
            } //UG
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[13][14];
            } //UU
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[13][15];
            } else {
                throw new RuntimeException("Invalid code for internal loop 1x1");
            }
        } /*GC-GC*/ else if (iNucleotide == NUCLEOTIDE_CODE.G && jNucleotide == NUCLEOTIDE_CODE.C && pNucleotide == NUCLEOTIDE_CODE.G && qNucleotide == NUCLEOTIDE_CODE.C) {
            //AA
            if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[14][0];
            } //AC
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[14][1];
            } //AG
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[14][2];
            } //AU
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[14][3];
            } //CA
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[14][4];
            } //CC
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[14][5];
            } //CG
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[14][6];
            } //CU
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[14][7];
            } //GA
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[14][8];
            } //GC
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[14][9];
            } //GG
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[14][10];
            } //GU
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[14][11];
            } //UA
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[14][12];
            } //UC
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[14][13];
            } //UG
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[14][14];
            } //UU
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[14][15];
            } else {
                throw new RuntimeException("Invalid code for internal loop 1x1");
            }
        } /*GC-UA*/ else if (iNucleotide == NUCLEOTIDE_CODE.G && jNucleotide == NUCLEOTIDE_CODE.C && pNucleotide == NUCLEOTIDE_CODE.U && qNucleotide == NUCLEOTIDE_CODE.A) {
            //AA
            if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[15][0];
            } //AC
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[15][1];
            } //AG
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[15][2];
            } //AU
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[15][3];
            } //CA
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[15][4];
            } //CC
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[15][5];
            } //CG
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[15][6];
            } //CU
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[15][7];
            } //GA
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[15][8];
            } //GC
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[15][9];
            } //GG
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[15][10];
            } //GU
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[15][11];
            } //UA
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[15][12];
            } //UC
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[15][13];
            } //UG
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[15][14];
            } //UU
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[15][15];
            } else {
                throw new RuntimeException("Invalid code for internal loop 1x1");
            }
        } /*GC-GU*/ else if (iNucleotide == NUCLEOTIDE_CODE.G && jNucleotide == NUCLEOTIDE_CODE.C && pNucleotide == NUCLEOTIDE_CODE.G && qNucleotide == NUCLEOTIDE_CODE.U) {
            //AA
            if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[16][0];
            } //AC
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[16][1];
            } //AG
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[16][2];
            } //AU
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[16][3];
            } //CA
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[16][4];
            } //CC
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[16][5];
            } //CG
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[16][6];
            } //CU
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[16][7];
            } //GA
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[16][8];
            } //GC
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[16][9];
            } //GG
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[16][10];
            } //GU
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[16][11];
            } //UA
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[16][12];
            } //UC
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[16][13];
            } //UG
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[16][14];
            } //UU
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[16][15];
            } else {
                throw new RuntimeException("Invalid code for internal loop 1x1");
            }
        } /*GC-UG*/ else if (iNucleotide == NUCLEOTIDE_CODE.G && jNucleotide == NUCLEOTIDE_CODE.C && pNucleotide == NUCLEOTIDE_CODE.U && qNucleotide == NUCLEOTIDE_CODE.G) {
            //AA
            if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[17][0];
            } //AC
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[17][1];
            } //AG
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[17][2];
            } //AU
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[17][3];
            } //CA
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[17][4];
            } //CC
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[17][5];
            } //CG
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[17][6];
            } //CU
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[17][7];
            } //GA
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[17][8];
            } //GC
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[17][9];
            } //GG
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[17][10];
            } //GU
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[17][11];
            } //UA
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[17][12];
            } //UC
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[17][13];
            } //UG
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[17][14];
            } //UU
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[17][15];
            } else {
                throw new RuntimeException("Invalid code for internal loop 1x1");
            }
        } //UA
        /*UA-AU*/ else if (iNucleotide == NUCLEOTIDE_CODE.U && jNucleotide == NUCLEOTIDE_CODE.A && pNucleotide == NUCLEOTIDE_CODE.A && qNucleotide == NUCLEOTIDE_CODE.U) {
            //AA
            if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[18][0];
            } //AC
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[18][1];
            } //AG
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[18][2];
            } //AU
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[18][3];
            } //CA
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[18][4];
            } //CC
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[18][5];
            } //CG
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[18][6];
            } //CU
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[18][7];
            } //GA
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[18][8];
            } //GC
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[18][9];
            } //GG
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[18][10];
            } //GU
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[18][11];
            } //UA
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[18][12];
            } //UC
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[18][13];
            } //UG
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[18][14];
            } //UU
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[18][15];
            } else {
                throw new RuntimeException("Invalid code for internal loop 1x1");
            }
        } /*UA-CG*/ else if (iNucleotide == NUCLEOTIDE_CODE.U && jNucleotide == NUCLEOTIDE_CODE.A && pNucleotide == NUCLEOTIDE_CODE.C && qNucleotide == NUCLEOTIDE_CODE.G) {
            //AA
            if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[19][0];
            } //AC
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[19][1];
            } //AG
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[19][2];
            } //AU
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[19][3];
            } //CA
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[19][4];
            } //CC
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[19][5];
            } //CG
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[19][6];
            } //CU
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[19][7];
            } //GA
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[19][8];
            } //GC
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[19][9];
            } //GG
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[19][10];
            } //GU
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[19][11];
            } //UA
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[19][12];
            } //UC
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[19][13];
            } //UG
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[19][14];
            } //UU
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[19][15];
            } else {
                throw new RuntimeException("Invalid code for internal loop 1x1");
            }
        } /*UA-GC*/ else if (iNucleotide == NUCLEOTIDE_CODE.U && jNucleotide == NUCLEOTIDE_CODE.A && pNucleotide == NUCLEOTIDE_CODE.G && qNucleotide == NUCLEOTIDE_CODE.C) {
            //AA
            if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[20][0];
            } //AC
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[20][1];
            } //AG
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[20][2];
            } //AU
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[20][3];
            } //CA
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[20][4];
            } //CC
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[20][5];
            } //CG
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[20][6];
            } //CU
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[20][7];
            } //GA
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[20][8];
            } //GC
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[20][9];
            } //GG
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[20][10];
            } //GU
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[20][11];
            } //UA
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[20][12];
            } //UC
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[20][13];
            } //UG
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[20][14];
            } //UU
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[20][15];
            } else {
                throw new RuntimeException("Invalid code for internal loop 1x1");
            }
        } /*UA-UA*/ else if (iNucleotide == NUCLEOTIDE_CODE.U && jNucleotide == NUCLEOTIDE_CODE.A && pNucleotide == NUCLEOTIDE_CODE.U && qNucleotide == NUCLEOTIDE_CODE.A) {
            //AA
            if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[21][0];
            } //AC
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[21][1];
            } //AG
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[21][2];
            } //AU
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[21][3];
            } //CA
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[21][4];
            } //CC
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[21][5];
            } //CG
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[21][6];
            } //CU
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[21][7];
            } //GA
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[21][8];
            } //GC
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[21][9];
            } //GG
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[21][10];
            } //GU
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[21][11];
            } //UA
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[21][12];
            } //UC
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[21][13];
            } //UG
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[21][14];
            } //UU
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[21][15];
            } else {
                throw new RuntimeException("Invalid code for internal loop 1x1");
            }
        } /*UA-GU*/ else if (iNucleotide == NUCLEOTIDE_CODE.U && jNucleotide == NUCLEOTIDE_CODE.A && pNucleotide == NUCLEOTIDE_CODE.G && qNucleotide == NUCLEOTIDE_CODE.U) {
            //AA
            if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[22][0];
            } //AC
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[22][1];
            } //AG
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[22][2];
            } //AU
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[22][3];
            } //CA
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[22][4];
            } //CC
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[22][5];
            } //CG
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[22][6];
            } //CU
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[22][7];
            } //GA
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[22][8];
            } //GC
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[22][9];
            } //GG
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[22][10];
            } //GU
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[22][11];
            } //UA
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[22][12];
            } //UC
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[22][13];
            } //UG
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[22][14];
            } //UU
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[22][15];
            } else {
                throw new RuntimeException("Invalid code for internal loop 1x1");
            }
        } /*UA-UG*/ else if (iNucleotide == NUCLEOTIDE_CODE.U && jNucleotide == NUCLEOTIDE_CODE.A && pNucleotide == NUCLEOTIDE_CODE.U && qNucleotide == NUCLEOTIDE_CODE.G) {
            //AA
            if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[23][0];
            } //AC
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[23][1];
            } //AG
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[23][2];
            } //AU
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[23][3];
            } //CA
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[23][4];
            } //CC
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[23][5];
            } //CG
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[23][6];
            } //CU
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[23][7];
            } //GA
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[23][8];
            } //GC
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[23][9];
            } //GG
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[23][10];
            } //GU
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[23][11];
            } //UA
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[23][12];
            } //UC
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[23][13];
            } //UG
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[23][14];
            } //UU
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[23][15];
            } else {
                throw new RuntimeException("Invalid code for internal loop 1x1");
            }
        } //GU
        /*GU-AU*/ else if (iNucleotide == NUCLEOTIDE_CODE.G && jNucleotide == NUCLEOTIDE_CODE.U && pNucleotide == NUCLEOTIDE_CODE.A && qNucleotide == NUCLEOTIDE_CODE.U) {
            //AA
            if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[24][0];
            } //AC
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[24][1];
            } //AG
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[24][2];
            } //AU
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[24][3];
            } //CA
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[24][4];
            } //CC
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[24][5];
            } //CG
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[24][6];
            } //CU
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[24][7];
            } //GA
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[24][8];
            } //GC
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[24][9];
            } //GG
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[24][10];
            } //GU
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[24][11];
            } //UA
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[24][12];
            } //UC
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[24][13];
            } //UG
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[24][14];
            } //UU
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[24][15];
            } else {
                throw new RuntimeException("Invalid code for internal loop 1x1");
            }
        } /*UA-CG*/ else if (iNucleotide == NUCLEOTIDE_CODE.G && jNucleotide == NUCLEOTIDE_CODE.U && pNucleotide == NUCLEOTIDE_CODE.C && qNucleotide == NUCLEOTIDE_CODE.G) {
            //AA
            if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[25][0];
            } //AC
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[25][1];
            } //AG
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[25][2];
            } //AU
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[25][3];
            } //CA
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[25][4];
            } //CC
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[25][5];
            } //CG
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[25][6];
            } //CU
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[25][7];
            } //GA
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[25][8];
            } //GC
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[25][9];
            } //GG
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[25][10];
            } //GU
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[25][11];
            } //UA
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[25][12];
            } //UC
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[25][13];
            } //UG
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[25][14];
            } //UU
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[25][15];
            } else {
                throw new RuntimeException("Invalid code for internal loop 1x1");
            }
        } /*UA-GC*/ else if (iNucleotide == NUCLEOTIDE_CODE.G && jNucleotide == NUCLEOTIDE_CODE.U && pNucleotide == NUCLEOTIDE_CODE.G && qNucleotide == NUCLEOTIDE_CODE.C) {
            //AA
            if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[26][0];
            } //AC
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[26][1];
            } //AG
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[26][2];
            } //AU
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[26][3];
            } //CA
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[26][4];
            } //CC
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[26][5];
            } //CG
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[26][6];
            } //CU
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[26][7];
            } //GA
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[26][8];
            } //GC
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[26][9];
            } //GG
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[26][10];
            } //GU
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[26][11];
            } //UA
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[26][12];
            } //UC
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[26][13];
            } //UG
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[26][14];
            } //UU
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[26][15];
            } else {
                throw new RuntimeException("Invalid code for internal loop 1x1");
            }
        } /*UA-UA*/ else if (iNucleotide == NUCLEOTIDE_CODE.G && jNucleotide == NUCLEOTIDE_CODE.U && pNucleotide == NUCLEOTIDE_CODE.U && qNucleotide == NUCLEOTIDE_CODE.A) {
            //AA
            if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[27][0];
            } //AC
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[27][1];
            } //AG
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[27][2];
            } //AU
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[27][3];
            } //CA
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[27][4];
            } //CC
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[27][5];
            } //CG
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[27][6];
            } //CU
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[27][7];
            } //GA
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[27][8];
            } //GC
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[27][9];
            } //GG
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[27][10];
            } //GU
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[27][11];
            } //UA
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[27][12];
            } //UC
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[27][13];
            } //UG
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[27][14];
            } //UU
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[27][15];
            } else {
                throw new RuntimeException("Invalid code for internal loop 1x1");
            }
        } /*UA-GU*/ else if (iNucleotide == NUCLEOTIDE_CODE.G && jNucleotide == NUCLEOTIDE_CODE.U && pNucleotide == NUCLEOTIDE_CODE.G && qNucleotide == NUCLEOTIDE_CODE.U) {
            //AA
            if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[28][0];
            } //AC
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[28][1];
            } //AG
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[28][2];
            } //AU
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[28][3];
            } //CA
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[28][4];
            } //CC
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[28][5];
            } //CG
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[28][6];
            } //CU
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[28][7];
            } //GA
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[28][8];
            } //GC
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[28][9];
            } //GG
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[28][10];
            } //GU
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[28][11];
            } //UA
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[28][12];
            } //UC
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[28][13];
            } //UG
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[28][14];
            } //UU
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[28][15];
            } else {
                throw new RuntimeException("Invalid code for internal loop 1x1");
            }
        } /*UA-UG*/ else if (iNucleotide == NUCLEOTIDE_CODE.G && jNucleotide == NUCLEOTIDE_CODE.U && pNucleotide == NUCLEOTIDE_CODE.U && qNucleotide == NUCLEOTIDE_CODE.G) {
            //AA
            if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[29][0];
            } //AC
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[29][1];
            } //AG
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[29][2];
            } //AU
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[29][3];
            } //CA
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[29][4];
            } //CC
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[29][5];
            } //CG
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[29][6];
            } //CU
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[29][7];
            } //GA
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[29][8];
            } //GC
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[29][9];
            } //GG
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[29][10];
            } //GU
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[29][11];
            } //UA
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[29][12];
            } //UC
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[29][13];
            } //UG
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[29][14];
            } //UU
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[29][15];
            } else {
                throw new RuntimeException("Invalid code for internal loop 1x1");
            }
        } //UG
        /*UG-AU*/ else if (iNucleotide == NUCLEOTIDE_CODE.U && jNucleotide == NUCLEOTIDE_CODE.G && pNucleotide == NUCLEOTIDE_CODE.A && qNucleotide == NUCLEOTIDE_CODE.U) {
            //AA
            if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[30][0];
            } //AC
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[30][1];
            } //AG
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[30][2];
            } //AU
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[30][3];
            } //CA
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[30][4];
            } //CC
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[30][5];
            } //CG
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[30][6];
            } //CU
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[30][7];
            } //GA
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[30][8];
            } //GC
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[30][9];
            } //GG
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[30][10];
            } //GU
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[30][11];
            } //UA
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[30][12];
            } //UC
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[30][13];
            } //UG
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[30][14];
            } //UU
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[30][15];
            } else {
                throw new RuntimeException("Invalid code for internal loop 1x1");
            }
        } /*UA-CG*/ else if (iNucleotide == NUCLEOTIDE_CODE.U && jNucleotide == NUCLEOTIDE_CODE.G && pNucleotide == NUCLEOTIDE_CODE.C && qNucleotide == NUCLEOTIDE_CODE.G) {
            //AA
            if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[31][0];
            } //AC
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[31][1];
            } //AG
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[31][2];
            } //AU
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[31][3];
            } //CA
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[31][4];
            } //CC
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[31][5];
            } //CG
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[31][6];
            } //CU
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[31][7];
            } //GA
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[31][8];
            } //GC
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[31][9];
            } //GG
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[31][10];
            } //GU
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[31][11];
            } //UA
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[31][12];
            } //UC
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[31][13];
            } //UG
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[31][14];
            } //UU
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[31][15];
            } else {
                throw new RuntimeException("Invalid code for internal loop 1x1");
            }
        } /*UA-GC*/ else if (iNucleotide == NUCLEOTIDE_CODE.U && jNucleotide == NUCLEOTIDE_CODE.G && pNucleotide == NUCLEOTIDE_CODE.G && qNucleotide == NUCLEOTIDE_CODE.C) {
            //AA
            if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[32][0];
            } //AC
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[32][1];
            } //AG
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[32][2];
            } //AU
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[32][3];
            } //CA
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[32][4];
            } //CC
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[32][5];
            } //CG
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[32][6];
            } //CU
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[32][7];
            } //GA
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[32][8];
            } //GC
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[32][9];
            } //GG
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[32][10];
            } //GU
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[32][11];
            } //UA
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[32][12];
            } //UC
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[32][13];
            } //UG
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[32][14];
            } //UU
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[32][15];
            } else {
                throw new RuntimeException("Invalid code for internal loop 1x1");
            }
        } /*UA-UA*/ else if (iNucleotide == NUCLEOTIDE_CODE.U && jNucleotide == NUCLEOTIDE_CODE.G && pNucleotide == NUCLEOTIDE_CODE.U && qNucleotide == NUCLEOTIDE_CODE.A) {
            //AA
            if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[33][0];
            } //AC
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[33][1];
            } //AG
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[33][2];
            } //AU
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[33][3];
            } //CA
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[33][4];
            } //CC
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[33][5];
            } //CG
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[33][6];
            } //CU
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[33][7];
            } //GA
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[33][8];
            } //GC
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[33][9];
            } //GG
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[33][10];
            } //GU
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[33][11];
            } //UA
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[33][12];
            } //UC
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[33][13];
            } //UG
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[33][14];
            } //UU
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[33][15];
            } else {
                throw new RuntimeException("Invalid code for internal loop 1x1");
            }
        } /*UA-GU*/ else if (iNucleotide == NUCLEOTIDE_CODE.U && jNucleotide == NUCLEOTIDE_CODE.G && pNucleotide == NUCLEOTIDE_CODE.G && qNucleotide == NUCLEOTIDE_CODE.U) {
            //AA
            if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[34][0];
            } //AC
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[34][1];
            } //AG
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[34][2];
            } //AU
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[34][3];
            } //CA
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[34][4];
            } //CC
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[34][5];
            } //CG
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[34][6];
            } //CU
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[34][7];
            } //GA
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[34][8];
            } //GC
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[34][9];
            } //GG
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[34][10];
            } //GU
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[34][11];
            } //UA
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[34][12];
            } //UC
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[34][13];
            } //UG
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[34][14];
            } //UU
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[34][15];
            } else {
                throw new RuntimeException("Invalid code for internal loop 1x1");
            }
        } /*UA-UG*/ else if (iNucleotide == NUCLEOTIDE_CODE.U && jNucleotide == NUCLEOTIDE_CODE.G && pNucleotide == NUCLEOTIDE_CODE.U && qNucleotide == NUCLEOTIDE_CODE.G) {
            //AA
            if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[35][0];
            } //AC
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[35][1];
            } //AG
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[35][2];
            } //AU
            else if (xNucleotide == NUCLEOTIDE_CODE.A && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[35][3];
            } //CA
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[35][4];
            } //CC
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[35][5];
            } //CG
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[35][6];
            } //CU
            else if (xNucleotide == NUCLEOTIDE_CODE.C && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[35][7];
            } //GA
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[35][8];
            } //GC
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[35][9];
            } //GG
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[35][10];
            } //GU
            else if (xNucleotide == NUCLEOTIDE_CODE.G && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[35][11];
            } //UA
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.A) {
                return ENERGY_CONSTANTS.internal_loop_11[35][12];
            } //UC
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.C) {
                return ENERGY_CONSTANTS.internal_loop_11[35][13];
            } //UG
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.G) {
                return ENERGY_CONSTANTS.internal_loop_11[35][14];
            } //UU
            else if (xNucleotide == NUCLEOTIDE_CODE.U && yNucleotide == NUCLEOTIDE_CODE.U) {
                return ENERGY_CONSTANTS.internal_loop_11[35][15];
            } else {
                throw new RuntimeException("Invalid code for internal loop 1x1");
            }
        } else {
            throw new RuntimeException("Invalid code for internal loop 1x1");
        }
    }
    
}
