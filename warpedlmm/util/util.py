import pysnptools.util.pheno 
from pysnptools.snpreader import Bed, Ped
import pysnptools.util as srutil
import numpy as np
import pandas


def get_snp_reader(snp_file):

    # get correct  SNP reader
    snp_file, ext_file = snp_file[:-4], snp_file[-4:]

    match ext_file:
        case '.bed':
            return Bed(snp_file)
        case '.ped':
            return Ped(snp_file)
        case _:
            raise ValueError('wrong extension for input SNP file')


def load_data(snp_file, pheno_file, covar_file, dist_file=None):

    # Load SNP data
    snp_reader = get_snp_reader(snp_file)


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

    # load distance file

    X = 1./np.sqrt((snp_data.val**2).sum() / float(snp_data.iid_count)) * snp_data.val

    if dist_file is None:

        K = np.dot(X, X.T) # TODO use symmetric dot to speed this up

    else:

        dist_reader = get_snp_reader(dist_file)
        dist_data = dist_reader.read().standardize()

        D = 1./np.sqrt((dist_data.val**2).sum() / float(dist_data.iid_count)) * dist_data.val
        K = np.dot(D, D.T) # TODO use symmetric dot to speed this up

    assert np.all(pheno['iid'] == snp_data.iid), "the samples are not sorted"

    return snp_data, pheno, covar, X, Y, K

def write_results_to_file(snp_data, pv, results_filename):
    results = pandas.DataFrame(index=snp_data.sid, columns=['Chr', 'ChrPos', 'Dist', 'PValue'])

    results['Chr'] = snp_data.pos[:, 0]
    results['Dist'] = snp_data.pos[:, 1]
    results['ChrPos'] = snp_data.pos[:, 2]
    results['PValue'] = pv[:, None]

    assert np.all(results.index == snp_data.sid) and np.all(results['PValue'] == pv), "the pvalues and/or SNP ids are not in order in the output file"

    results.sort_values(by='PValue', inplace=True)
    results['warped'] = 1
    results.to_feather(results_filename)
