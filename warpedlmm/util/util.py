import pysnptools.util.pheno 
from pysnptools.snpreader import Bed, Ped
import pysnptools.util as srutil
import numpy as np
import pandas

def load_data(snp_file, pheno_file, covar_file):
    # Load SNP data
    snp_file, ext_file = snp_file[:-4], snp_file[-4:]

    match ext_file:
        case '.bed':
            snp_reader = Bed(snp_file)
        case '.ped':
            snp_reader = Ped(snp_file)
        case _:
            raise ValueError('wrong extension for input SNP file')

    # Load phenotype
    pheno = pysnptools.util.pheno.loadPhen(pheno_file)

    # Load covariates
    if covar_file is not None:
        covar = pysnptools.util.pheno.loadPhen(covar_file)
        snp_reader, pheno, covar = srutil.intersect_apply([snp_reader, pheno, covar])
        covar = covar['vals']
    else:
        snp_reader, pheno = srutil.intersect_apply([snp_reader, pheno])
        covar = None

    snp_data = snp_reader.read().standardize()
    Y = pheno['vals']
    Y -= Y.mean(0)
    Y /= Y.std(0)

    X = 1./np.sqrt((snp_data.val**2).sum() / float(snp_data.iid_count)) * snp_data.val
    K = np.dot(X, X.T) # TODO use symmetric dot to speed this up

    assert np.all(pheno['iid'] == snp_data.iid), "the samples are not sorted"

    return snp_data, pheno, covar, X, Y, K

def write_results_to_file(snp_data, pv, results_filename):
    results = pandas.DataFrame(index=snp_data.sid, columns=['Chr', 'ChrPos', 'Dist', 'PValue'])

    results['Chr'] = snp_data.pos[:, 0]
    results['Dist'] = snp_data.pos[:, 1]
    results['ChrPos'] = snp_data.pos[:, 2]
    results['PValue'] = pv[:, None]

    assert np.all(results.index == snp_data.sid) and np.all(results['PValue'] == pv), "the pvalues and/or SNP ids are not in order in the output file"

    results.to_feather(results_filename)
