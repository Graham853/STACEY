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

import beast.core.Citation;
import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.speciation.SpeciesTreeDistribution;
import beast.evolution.tree.Node;
import beast.evolution.tree.TreeInterface;

import java.util.List;



@Description("Birth-death-collapse distribution of prior for species tree or minimal clusters tree.")

@Citation(value = "Jones and Oxelman. " +
        "DISSECT: an assignment-free Bayesian discovery method for species delimitation under the multispecies coalescent. " +
        "http://biorxiv.org/content/early/2014/03/03/003178. doi: 10.1101/003178.",
        firstAuthorSurname = "Jones")


public class BirthDeathCollapseModel extends SpeciesTreeDistribution {
    @SuppressWarnings("WeakerAccess")
    public static final Input<Double> collapseHeight =
            new Input<Double>("collapseHeight",
                    "Collapse height value, epsilon in birth/death/collapse model.",
                    Input.Validate.REQUIRED);

    @SuppressWarnings({"CanBeFinal", "WeakerAccess"})
    public Input<RealParameter> birthDiffRate =
            new Input<RealParameter>("birthDiffRate",
                    "Growth rate rate parameter, lambda-mu in birth/death/collapse model.",
                    Input.Validate.REQUIRED);

    @SuppressWarnings({"CanBeFinal", "WeakerAccess"})
    public Input<RealParameter> relativeDeathRate =
            new Input<RealParameter>("relativeDeathRate",
                    "Relative death rate parameter, mu/lambda in birth/death/collapse model.",
                    Input.Validate.REQUIRED);

    @SuppressWarnings({"CanBeFinal", "WeakerAccess"})
    public Input<RealParameter> collapseWeight =
            new Input<RealParameter>("collapseWeight",
                    "Collapse weight parameter, w in birth/death/collapse model.",
                    Input.Validate.REQUIRED);

    @SuppressWarnings({"CanBeFinal", "WeakerAccess"})
    public Input<RealParameter> originHeight =
            new Input<RealParameter>("originHeight",
                    "Origin height of the species-or-minimal-clusters tree. Must be estimated. Initial value is ignored. In birth/death/collapse model",
                    Input.Validate.REQUIRED);



    @Override
    public void initAndValidate() throws Exception {
        super.initAndValidate();

        if ((collapseHeight.get() < 1e-30) || (collapseHeight.get() > 1e10)) {
            throw new Exception("Bad collapseHeight value");
        }
        /*//TODO grjtodo. This fails if w is fixed. No lower/upper value?
        if (collapseWeight.get().lowerValueInput.get() < 0  ||
                collapseWeight.get().upperValueInput.get() > 1.0 ) {
            throw new Exception("Bad collapseWeight limits");
        }
        if (relativeDeathRate.get().lowerValueInput.get() < 0  ||
                relativeDeathRate.get().upperValueInput.get() > 1.0 ) {
            throw new Exception("Bad relativeDeathRate limits");
        }*/
        // Initialise originHeight to be a bit bigger than root height
        TreeInterface tree = treeInput.get();
        Double [] initOHarray = new Double[1];
        initOHarray[0] = 1.05 * tree.getRoot().getHeight();
        RealParameter initOH = new RealParameter(initOHarray);
        originHeight.get().assignFromWithoutID(initOH);
    }



    @Override
    public double calculateTreeLogLikelihood(final TreeInterface tree) {
        logP = treeLL(tree);
        return logP;
    }



    private double treeLL(final TreeInterface tree) {
        double logpt = 0.0;
        int ntips = tree.getLeafNodeCount();
        double alpha = birthDiffRate.get().getValue();
        double beta = relativeDeathRate.get().getValue();
        double tor = originHeight.get().getValue();
        double w = collapseWeight.get().getValue();

        double rooth = tree.getRoot().getHeight();
        if (rooth > tor) {
            return Double.NEGATIVE_INFINITY;
        }

        logpt += originHeightLogLikelihood(tor, alpha, beta, w, ntips);

        List<Node> internalNodes = tree.getInternalNodes();
        for (Node node : internalNodes) {
            final double height = node.getHeight();
            double usualpn = nodeHeightLikelihood(height, tor, alpha, beta);
            double collapsedpn = belowCollapseHeight(height) ? (1.0 / collapseHeight.get()) : 0.0;
            logpt += Math.log((1.0 - w) * usualpn + w * collapsedpn);
        }
        return logpt;
    }


    // provided to help avoid inconsistent treatment of h == collapseHeight
    static boolean belowCollapseHeight(double h) {
        return (h < collapseHeight.get());
    }


   static public double collapseHeight() {
        return collapseHeight.get().doubleValue();
    }


    private double originHeightLogLikelihood(double t, double a, double b, double w, int n) {
        double E = Math.exp(-a * t);
        double B = (1 - E) / (1-b*E);
        double z = 0.0;
        z += Math.log(a);
        z += Math.log(1 - b);
        z -= a * t;
        z -= 2 * Math.log(1 - b * E);
        z += (n-2) * Math.log(w + (1 - w) * B);
        z +=  Math.log(w + n * (1 - w) * B);
        return z;
    }


    private double nodeHeightLikelihood(double s, double t, double a, double b) {
        double Es = Math.exp(-a * s);
        double Et = Math.exp(-a * t);
        double z = 0.0;
        if (s < t) {
            z = a;
            z *= (1 - b);
            z *= Es;
            z /= (1 - b * Es) * (1 - b * Es);
            z *= (1 - b * Et);
            z /= (1 - Et);
        }
        return z;
    }

    @Override
    protected boolean requiresRecalculation() {
        return super.requiresRecalculation() ||
                birthDiffRate.get().somethingIsDirty() ||
                relativeDeathRate.get().somethingIsDirty() ||
                collapseWeight.get().somethingIsDirty() ||
                originHeight.get().somethingIsDirty();
    }
}
