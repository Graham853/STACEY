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

import java.util.*;


@Description("Imporatant utility class for STACEY, used by several operators. " +
        "It adds a set of species (or minimal clusters) to every node in the SMC-tree and in all gene trees.")
public class UnionArrays {

    private BitUnion [] sUnions;
    private BitUnion [][] gUnions;

    private TreeInterface sTree;
    private List<Tree> gTrees;
    private static UnionArrays unionArrays = null;


    public static UnionArrays initialise(TreeInterface sTree, List<Tree> gTrees, Bindings bindings) {
        if (unionArrays == null) {
            unionArrays = new UnionArrays(sTree, gTrees, bindings);
        }
        return unionArrays;
    }



    public void update() {
        beastSubtreeToUnions(sUnions, sTree.getRoot());
        for (int j = 0; j < gUnions.length; j++) {
            beastSubtreeToUnions(gUnions[j], gTrees.get(j).getRoot());
        }
    }



    public void reset() {
        // needed?
    }


    public BitUnion sNodeUnion(int n) {
        return sUnions[n];
    }


    public BitUnion gNodeUnion(int j, int i) {
        return gUnions[j][i];
    }


    public ArrayList<Node> getStraddlers(int j, BitUnion x, BitUnion y) {
        ArrayList<Node> s = new ArrayList<>();
        subtreeStraddlers(s, j, x, y, gTrees.get(j).getRoot());
        return s;
    }



    // for SPR move, maybe can avoid there
    public int nodeIndexOfUnionInSubSTree(Node node, BitUnion x) {
        if (node.isLeaf()) {
            return node.getNr();
        }
        Node lftNode = node.getChild(0);
        Node rgtNode = node.getChild(1);
        if (x.isContainedIn(sUnions[lftNode.getNr()])) {
            return nodeIndexOfUnionInSubSTree(lftNode, x);
        } else if (x.isContainedIn(sUnions[rgtNode.getNr()])) {
            return nodeIndexOfUnionInSubSTree(rgtNode, x);
        } else {
            return node.getNr();
        }
    }


    /*******************************************************************************/


    private UnionArrays(TreeInterface sTree, List<Tree> gTrees, Bindings bindings) {
        this.sTree = sTree;
        this.gTrees = gTrees;

        sUnions = new BitUnion[sTree.getNodeCount()];
        int nSMCs = bindings.smcTreeTipCount();
        for (int n = 0; n < sUnions.length; n++) {
            sUnions[n] = new BitUnion(nSMCs);
            if (n < nSMCs) {
                sUnions[n].replaceWith(bindings.tipUnionOfSMCNodeNr(n));
            }
        }
        gUnions = new BitUnion[gTrees.size()][];
        for (int j = 0; j < gUnions.length; j++) {
            gUnions[j] = new BitUnion[gTrees.get(j).getNodeCount()];
            int nGTips = gTrees.get(j).getLeafNodeCount();
            for (int i = 0; i < gUnions[j].length; i++) {
                gUnions[j][i] = new BitUnion(nSMCs);
                if (i < nGTips) {
                    gUnions[j][i].replaceWith(bindings.tipUnionOfGNodeNr(j, i));
                }
            }
        }
    }



    private void beastSubtreeToUnions(BitUnion [] unions, Node node) {
        Node lftNode = node.getChild(0);
        if (!lftNode.isLeaf()) {
            beastSubtreeToUnions(unions, lftNode);
        }
        Node rgtNode = node.getChild(1);
        if (!rgtNode.isLeaf()) {
            beastSubtreeToUnions(unions, rgtNode);
        }
        unions[node.getNr()].replaceWith(unions[lftNode.getNr()]);
        unions[node.getNr()].union(unions[rgtNode.getNr()]);
    }



    private void subtreeStraddlers(ArrayList<Node> s, int j, BitUnion x, BitUnion y, Node node) {
        Node lftN = node.getChild(0);
        Node rgtN = node.getChild(1);
        boolean lftx = gUnions[j][lftN.getNr()].overlaps(x);
        boolean lfty = gUnions[j][lftN.getNr()].overlaps(y);
        boolean rgtx = gUnions[j][rgtN.getNr()].overlaps(x);
        boolean rgty = gUnions[j][rgtN.getNr()].overlaps(y);
        if ((lftx & rgty) | (lfty & rgtx)) {
            s.add(node);
        }
        if (lftx & lfty) {
            subtreeStraddlers(s, j, x, y, lftN);
        }
        if (rgtx & rgty) {
            subtreeStraddlers(s, j, x, y, rgtN);
        }
    }

}
