/*
 * Copyright (c) 2010 The Jackson Laboratory
 * 
 * This is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this software.  If not, see <http://www.gnu.org/licenses/>.
 */
package org.jax.emma;

/**
 * for filtering data by sex
 */
public enum SexFilter
{
    /**
     * allow females through
     */
    ALLOW_FEMALE
    {
        /**
         * {@inheritDoc}
         */
        @Override
        public String toString()
        {
            return "Only Allow Females";
        }
    },
    
    /**
     * allow males through
     */
    ALLOW_MALE
    {
        /**
         * {@inheritDoc}
         */
        @Override
        public String toString()
        {
            return "Only Allow Males";
        }
    },
    
    /**
     * allow either female or male
     */
    AGNOSTIC
    {
        /**
         * {@inheritDoc}
         */
        @Override
        public String toString()
        {
            return "Sex Agnostic";
        }
    }
}