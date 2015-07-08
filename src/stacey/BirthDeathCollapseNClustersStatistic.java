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

package stacey;

import beast.core.BEASTObject;
import beast.core.Description;
import beast.core.Input;
import beast.core.Loggable;
import beast.evolution.tree.Node;
import beast.evolution.tree.TreeInterface;

import java.io.PrintStream;
import java.util.List;

/**
 *  Created by Graham Jones on 06/09/2014.
 */

@Description("Statistic that counts the number of clusters (equals one more than the number of collapsed nodes)")
public class BirthDeathCollapseNClustersStatistic  extends BEASTObject implements Loggable {
    @SuppressWarnings({"CanBeFinal", "WeakerAccess"})
    public Input<TreeInterface> smcTree =
            new Input<TreeInterface>("smcTree",
                    "Species or minimal clusters  tree", Input.Validate.REQUIRED);


    // bdcm not used here but it is in templates and XML files
    @SuppressWarnings("UnusedDeclaration")
    public Input<BirthDeathCollapseModel> bdcm =
            new Input<BirthDeathCollapseModel>("bdcm", "The birth death collapse model", Input.Validate.REQUIRED);



    @Override
    public void initAndValidate() {
    }



    @Override
    public void init(PrintStream out) throws Exception {
        out.print("NClusters\t");
    }

    @Override
    public void log(int nSample, PrintStream out) {
        out.print(""+getStatisticValue() + "\t");
    }

    @Override
    public void close(PrintStream out) {
        // nothing to do
    }

    private double getStatisticValue() {
        List<Node> inodes = smcTree.get().getInternalNodes();
        int n =  0;
        for (Node node : inodes) {
            double h = node.getHeight();
            if (!BirthDeathCollapseModel.belowCollapseHeight(h)) {
                n++;
            }
        }
        return n+1;
    }


}
