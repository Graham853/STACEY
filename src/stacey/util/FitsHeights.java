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
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.evolution.tree.TreeInterface;

import java.util.ArrayList;
import java.util.List;




@Description("Important utility class for STACEY, used for logP calculation in PIOMSCoalescentDistribution. " +
        "It connects gene tree nodes to SMC-tree nodes and branches.")
public class FitsHeights {

    private ArrayList<ArrayList<ArrayList<Double>>> heights;
    // heights.get(n).get(j) is an array of the coalescence heights from gtree j
    // which are inside stree branch n

    private Bindings bindings;

    private int currentGTreeIndex;
    private TreeInterface sTree;
    private List<Tree> gTrees;
    // I assume it is more efficient to store these once than to pass them down a recursive routine.



    /*************************************************************************/

    public FitsHeights(TreeInterface sTree, List<Tree> gTrees, Bindings bindings) {
        this.bindings = bindings;
        this.sTree = sTree;
        this.gTrees = gTrees;
        initHeightArray(sTree, gTrees);
    }


    // Returns false if gtree j doesn't fit, else updates a 3D array of heights.
    public boolean updateFitHeightsForOneGTree(int j) {
        currentGTreeIndex = j;
        // clear old heights
        for (ArrayList<ArrayList<Double>> hostNodeHeights: heights) {
            hostNodeHeights.get(j).clear();
        }
        Node node = makeHeightsSubtree(gTrees.get(j).getRoot());
        return (node != null);
    }

    // Once the heights array is updated, returns a 1D array of the coalescence
    // heights from gtree j which are inside stree branch n.
    public ArrayList<Double> getHeightsFromSNodeNrGTree(int n, int j) {
        return heights.get(n).get(j);
    }



    /************************  private  *******************************************/


    private void initHeightArray(TreeInterface sTree, List<Tree> gTrees) {
        int numSNodes = sTree.getNodeCount();
        int numGTrees = gTrees.size();
        heights = new ArrayList<>(numSNodes);
        for (int n = 0; n < numSNodes; n++) {
            heights.add(new ArrayList<>(numGTrees));
            for (int j = 0; j < numGTrees; j++) {
                heights.get(n).add(new ArrayList<>(0));
            }
        }
    }



    // This is speed-critical, eg 19 indiv in 9 spp, 100 loci, this took 13% of run time
    private Node makeHeightsSubtree(Node hereG) {
        Node lftG = hereG.getChild(0);
        Node lftS;
        if (lftG.isLeaf()) {
            lftS = sTree.getNode(bindings.smcTipNrFromGTreeTipNr(currentGTreeIndex, lftG.getNr()));
        } else {
            lftS = makeHeightsSubtree(lftG);
            if (lftS == null) {
                return null;
            }
        }
        Node rgtG = hereG.getChild(1);
        Node rgtS;
        if (rgtG.isLeaf()) {
            rgtS = sTree.getNode(bindings.smcTipNrFromGTreeTipNr(currentGTreeIndex, rgtG.getNr()));
        } else {
            rgtS = makeHeightsSubtree(rgtG);
            if (rgtS == null) {
                return null;
            }
        }
        while (lftS != rgtS) {
            if (lftS.getHeight() < rgtS.getHeight()) {
                lftS = lftS.getParent();
            } else {
                rgtS = rgtS.getParent();
            }
        }
        Node hereS = lftS;
        double heightG = hereG.getHeight();
        if (hereS.getHeight() > heightG) {
            return null;
        }
        Node hostS = hereS;
        while (!hostS.isRoot()  &&  hostS.getParent().getHeight() < heightG) {
            hostS = hostS.getParent();
        }
        heights.get(hostS.getNr()).get(currentGTreeIndex).add(heightG);
        return hereS;
    }
}
