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
import beast.core.parameter.RealParameter;

/**
 *  Created by Graham Jones on 20/08/2014.
 */

@Description("A component in a mixture of inverse gamma distributions.")


public class InverseGammaComponent extends BEASTObject {
    @SuppressWarnings({"CanBeFinal", "WeakerAccess"})
    public Input<RealParameter> weight =
            new Input<RealParameter>("weight",
                    "weight of component in mixture of inverse gammas", Input.Validate.REQUIRED);

    @SuppressWarnings({"CanBeFinal", "WeakerAccess"})
    public Input<RealParameter>  alpha =
            new Input<RealParameter>("alpha",
                    "shape parameter of an inverse gamma", Input.Validate.REQUIRED);

    @SuppressWarnings({"CanBeFinal", "WeakerAccess"})
    public Input<RealParameter> beta =
            new Input<RealParameter>("beta",
                    "scale parameter of an inverse gamma", Input.Validate.REQUIRED);

    // inv gamma pdf is parameterized as  b^a/Gamma(a)  x^(-a-1)  exp(-b/x)
    // mean is b/(a-1) if a>1, var is  b^2/((a-1)^2 (a-2)) if a>2.



    @Override
    public void initAndValidate() throws Exception {
        if (getWeight() <= 0.0) {
            throw new Exception("weight must be positive");
        }
        if (getAlpha() <= 0.0) {
            throw new Exception("alpha must be positive");
        }
        if (getBeta() <= 0.0) {
            throw new Exception("beta must be positive");
        }
    }

    double getWeight() { return weight.get().getValue(); }
    double getAlpha() { return alpha.get().getValue(); }
    double getBeta() { return beta.get().getValue(); }

    void normalizeWeight(double total) throws Exception {
        Double [] normWt = new Double[1];
        normWt[0] = weight.get().getValue() / total;
        RealParameter normalizedWeight = new RealParameter(normWt);
        weight.get().assignFromWithoutID(normalizedWeight);
    }

}
