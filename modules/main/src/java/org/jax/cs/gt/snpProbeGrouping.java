/*
* Copyright (c) 2011 The Jackson Laboratory
*
* This is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* This software is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this software. If not, see <http://www.gnu.org/licenses/>.
*/

package org.jax.cs.gt;

import java.util.ArrayList;
import java.util.List;

/**
 *
 * @author gbeane
 */
public class snpProbeGrouping {

    private List<Float> intensitiesA;
    private List<Float> intensitiesB;
    private String id;

    public snpProbeGrouping(String id)
    {
        this.id = id;
        intensitiesA = new ArrayList<Float>();
        intensitiesB = new ArrayList<Float>();
    }

    public void insertAIntensity(float intensity)
    {
        intensitiesA.add(intensity);
    }

    public void insertBIntensity(float intensity)
    {
        intensitiesB.add(intensity);
    }

    public List getAIntensities()
    {
        return intensitiesA;
    }

    public List getBIntensities()
    {
        return intensitiesB;
    }

    public String getID()
    {
        return id;
    }

    public static float mean(List<Float> intensities)
    {
        float sum = 0;
        for (Float val : intensities) {
            sum += val;
        }

        return sum / intensities.size();
    }

}
