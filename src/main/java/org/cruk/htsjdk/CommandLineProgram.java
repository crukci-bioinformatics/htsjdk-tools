/**
 * Copyright (c) 2021 CRUK Cambridge Institute - Bioinformatics Core
 *
 * Licensed under the MIT license (http://opensource.org/licenses/MIT)
 * This file may not be copied, modified, or distributed except according
 * to those terms.
 */

package org.cruk.htsjdk;

import picocli.CommandLine.IVersionProvider;

/**
 * Abstract base class for command line programs with functions obtaining
 * package and version information.
 *
 * @author eldrid01
 */
public abstract class CommandLineProgram implements Runnable, IVersionProvider {

    /**
     * @return the package name and version
     */
    protected String getPackageNameAndVersion() {
        return getClass().getPackage().getImplementationTitle() + " "
                + getClass().getPackage().getImplementationVersion();
    }

    @Override
    public String[] getVersion() throws Exception {
        return new String[] { getPackageNameAndVersion() + "  " + getClass().getName() };
    }
}