/*
        Copyright (C) 2015 Graham Jones, www.indriid.com

        This file is part of STACEY.

        STACEY is free software: you can redistribute it and/or modify
        it under the terms of the GNU General Public License as published by
        the Free Software Foundation, either version 3 of the License, or
        (at your option) any later version.

        STACEY is distributed in the hope that it will be useful,
        but WITHOUT ANY WARRANTY; without even the implied warranty of
        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
        GNU General Public License for more details.

        You should have received a copy of the GNU General Public License
        along with STACEY.  If not, see <http://www.gnu.org/licenses/>.
*/

package stacey.util;

import beast.core.Description;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.evolution.tree.TreeInterface;

import java.util.ArrayList;
import java.util.List;
import java.util.Stack;

/**
 *  Created by Graham Jones on 27/08/2014.
 */



/* Major rewrite, June 2015.
A Bindings is constructed initially from the smcTree and list of gTrees. It only stores
information which is constant for the analysis, mainly about the tips of the trees.

It serves two main purposes, for logP calculation of PIOMSCoalescentDistribution, and for operators.

For logP. It provides nLineagesForBeastTipNrAndGtree(n,j) which returns the number of lineages at gtree tip n
in gtree j, where n is a beast node nr for the tip.

For operators. It stores information about the tips of smcTree and gTrees, mapping indices of tips
of gtrees to BitUnions of the smcTree.
These can be used to make/update union arrays:
smcTreeTipCount()
tipUnionOfSMCNodeNr(int nodeNr)
tipUnionOfGNodeNr(int j, int nodeNr)
*/

@Description("Provides various linkages between the SMC-Tree tips and gene tree tips.")
public class Bindings {
    private final BitUnion [] smcTipUnions;
    private final BitUnion [][] gTipNrToUnion;
    // gTipNrToUnion[j][i] is the union for (tip) node i in gtree j,
    // where 0 <= i < number of tips in gtree j.
    private final int [][] gTipNrToSmcTipNr;
    // gTipNrToSmcTipNr[j][i] is the stree node number for (tip) node i in gtree j,
    // where 0 <= i < number of tips in gtree j. Can use this together with
    // stree to get the stree tip node to which the gtree node belongs.
    private final int [][] nLinsForSmcTipNr;
    // nLinsForSmcTipNr[n][j] is the number of lineages from gene tree j
    // which are assigned to the stree (tip) node n
    // where 0 <= n < number of tips in stree j.

    private static Bindings bindings = null;


    /************************* all clients call this before other methods ******************/

    public static Bindings initialise(TreeInterface sTree, List<Tree> gTrees) {
        if (bindings == null) {
            bindings = new Bindings(sTree, gTrees);
        }
        return bindings;
    }


    /***********************************  For logP  *********************************************************/


    // Used by FitsHeights
    public int [][] getGTipNrToSmcTipNr() {
        return gTipNrToSmcTipNr;
    }



    public int nLineagesForBeastTipNrAndGtree(int beastNodeNr, int j) {
        if (beastNodeNr >= nLinsForSmcTipNr.length) {
            System.out.println("DEBUGGING\n");
        }
        if (j >= nLinsForSmcTipNr[beastNodeNr].length) {
            System.out.println("DEBUGGING\n");
        }
        return nLinsForSmcTipNr[beastNodeNr][j];
    }



    /***********************************  For operators  **********************************************/

    public int smcTreeTipCount() {
        return smcTipUnions.length;
    }


    public BitUnion tipUnionOfSMCNodeNr(int nodeNr) {
        return smcTipUnions[nodeNr];
    }

    public BitUnion tipUnionOfGNodeNr(int j, int nodeNr) {
        return gTipNrToUnion[j][nodeNr];
    }


    /***********************************************************************************/
    /********************************* private *****************************************/
    /***********************************************************************************/



    private Bindings(TreeInterface sTree, List<Tree> gTrees) {
        int nSTips = sTree.getLeafNodeCount();
        int nGTrees = gTrees.size();
        smcTipUnions = new BitUnion[nSTips];
        // make tip unions for smc tree
        for (int i = 0; i < nSTips; i++) {
            smcTipUnions[i] = new BitUnion(nSTips);
            smcTipUnions[i].insert(i);
        }
        gTipNrToUnion = new BitUnion[nGTrees][];
        gTipNrToSmcTipNr = new int[nGTrees][];
        nLinsForSmcTipNr = new int[nSTips][nGTrees];
        setUpTipMaps(sTree, gTrees);
        setUpTipNlineages(sTree, gTrees);
    }

    private void setUpTipNlineages(TreeInterface sTree, List<Tree> gTrees) {
        int nSTips = sTree.getLeafNodeCount();
        int nGTrees = gTrees.size();
        for (int STipNr = 0; STipNr < nSTips; STipNr++) {
            for (int j = 0; j < nGTrees; j++) {
                int nlins = 0;
                for (int GTipNr = 0; GTipNr < gTipNrToUnion[j].length; GTipNr++) {
                    if (gTipNrToSmcTipNr[j][GTipNr] == STipNr) {
                        nlins++;
                    }
                }
                nLinsForSmcTipNr[STipNr][j] = nlins;
            }
        }
    }



    private void setUpTipMaps(TreeInterface sTree, List<Tree> gTrees) {
        for (int j = 0; j < gTrees.size(); j++) {
            TreeInterface gTree = gTrees.get(j);
            List<Node> gTips = gTree.getExternalNodes();
            int nGTips = gTips.size();
            gTipNrToUnion[j] = new BitUnion[nGTips];
            gTipNrToSmcTipNr[j] = new int[nGTips];
            for (int i = 0; i < nGTips; i++) {
                int beastNr = gTips.get(i).getNr();
                assert i == beastNr;
                String gtName = gTips.get(i).getID();
                int smcTipNr = smcTipNrFromGTipID(sTree, gtName);
                gTipNrToSmcTipNr[j][beastNr] = smcTipNr;
                gTipNrToUnion[j][beastNr] = new BitUnion(sTree.getLeafNodeCount());
                gTipNrToUnion[j][beastNr].insert(smcTipNr);
            }
        }
    }



    private int smcTipNrFromGTipID(TreeInterface sTree, String gTipID) {
        String smcTipID = smcTipIDFromGTipID(sTree.getTaxonset(), gTipID);
        List<Node> sTips = sTree.getExternalNodes();
        int tipNr = -1;
        for (Node sTip : sTips) {
            if (sTip.getID().compareTo(smcTipID) == 0) {
                if (tipNr != -1) {
                    System.out.println("Error in smcTipNrFromGTipID() with taxon '" + gTipID + "'.\n" +
                            "Assigned twice to SMC-tree tip number " + tipNr + ".");
                    System.exit(1);
                }
                assert tipNr == -1;
                tipNr = sTip.getNr();
            }
        }
        if (tipNr < 0) {
            System.out.println("Error in smcTipNrFromGTipID() with taxon '" + gTipID + "'.\n" +
            "Can't assign to SMC-tree tip number " + tipNr + ".");
            System.exit(1);
        }
        assert tipNr >= 0;
        return tipNr;
    }


    private String smcTipIDFromGTipID(TaxonSet taxonSetOfSets, String gTipID) {
        String smcTipID = "";
        for (int sts = 0; sts < taxonSetOfSets.getTaxonCount(); sts++) {
            TaxonSet smc = (TaxonSet) taxonSetOfSets.taxonsetInput.get().get(sts);
            for (int gt = 0; gt < smc.getTaxonCount(); gt++) {
                String sgID = smc.taxonsetInput.get().get(gt).getID();
                if (gTipID.compareTo(sgID) == 0) {
                    if (smcTipID.compareTo("") != 0) {
                        System.out.println("Error in smcTipIDFromGTipID() with taxon '" + gTipID + "'.\n" +
                                "Assigned twice to SMC-tree tip '" + smcTipID + "'.");
                        System.exit(1);
                    }
                    assert smcTipID.compareTo("") == 0;
                    smcTipID = smc.getID();
                }
            }
        }
        if (smcTipID.compareTo("") == 0) {
            System.out.println("Error in smcTipIDFromGTipID() with taxon '" + gTipID + "'.\n" +
                    "Can't assign to SMC-tree tip '" + smcTipID + "'.");
            System.exit(1);
        }
        assert smcTipID.compareTo("") != 0;
        return smcTipID;
    }

}
