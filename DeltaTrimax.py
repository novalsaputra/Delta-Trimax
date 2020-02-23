# -*- coding: utf-8 -*-

import numpy as np

class DeltaTrimax4():

    def __init__(self, D):
        self.D = D.copy()

    def _compute_MSR(self, time, gene, sample, gene_add = False, sample_add = False, time_add = False):
 
        time_idx = np.expand_dims(np.expand_dims(np.nonzero(time)[0], axis=1), axis=1)
        gene_idx = np.expand_dims(np.expand_dims(np.nonzero(gene)[0], axis=0), axis=2)
        sample_idx = np.expand_dims(np.expand_dims(np.nonzero(sample)[0], axis=0), axis=0)

        if (not time_idx.size) or (not gene_idx.size) or (not sample_idx.size):
            raise EmptyTriclusterException()

        subarr = self.D[time_idx, gene_idx, sample_idx]
        self.n_time = subarr.shape[0]
        self.n_gene = subarr.shape[1]
        self.n_sample = subarr.shape[2]


        # Computation of m_iJK
        if gene_add == True:
            pass
        else:
            m_iJK = np.nanmean(np.nanmean(subarr,axis=2),axis=0)
            self.m_iJK = np.expand_dims(np.expand_dims(m_iJK,axis=0),axis=2)

        # Computation of m_IjK
        if sample_add == True:
            pass
        else:
            m_IjK = np.nanmean(np.nanmean(subarr, axis=0), axis=0)
            self.m_IjK = np.expand_dims(np.expand_dims(m_IjK, axis=0), axis=0)

        # Computation of m_IJk
        if time_add == True:
            pass
        else:
            m_IJk = np.nanmean(np.nanmean(subarr,axis=1),axis=1)
            self.m_IJk = np.expand_dims(np.expand_dims(m_IJk,axis=1),axis=2)

        # Computation of m_IJK
        m_IJK = np.nanmean(subarr)

        # Computation of MSR
        residue = subarr - self.m_iJK - self.m_IjK - self.m_IJk + (2*m_IJK)
        SR = np.square(residue)

        self.MSR = np.nanmean(SR)
        self.MSR_time = np.nanmean(np.nanmean(SR, axis=2), axis=1)
        self.MSR_gene = np.nanmean(np.nanmean(SR, axis=2), axis=0)
        self.MSR_sample = np.nanmean(np.nanmean(SR, axis=0), axis=0)

        # Check tolerance
        self.MSR_time[self.MSR_time < self.tol] = 0
        self.MSR_gene[self.MSR_gene < self.tol] = 0
        self.MSR_sample[self.MSR_sample < self.tol] = 0
        self.MSR = 0 if (self.MSR < self.tol or np.isnan(self.MSR)) else self.MSR

    def _single_node_deletion(self, time, gene, sample):

        self._compute_MSR(time, gene, sample)

        while (self.MSR > self.delta):
            gene_idx = np.nanargmax(self.MSR_gene)
            sample_idx = np.nanargmax(self.MSR_sample)
            time_idx = np.nanargmax(self.MSR_time)

            if (self.MSR_gene[gene_idx] > self.MSR_sample[sample_idx]):
                if (self.MSR_gene[gene_idx] >  self.MSR_time[time_idx]):
                    # Delete Gene
                    nonz_idx = gene.nonzero()[0]
                    gene.put(nonz_idx[gene_idx],0)
                else:
                    # Delete Time
                    nonz_idx = time.nonzero()[0]
                    time.put(nonz_idx[time_idx],0)
            else:
                if (self.MSR_sample[sample_idx] > self.MSR_time[time_idx]):
                    # Delete Sample
                    nonz_idx = sample.nonzero()[0]
                    sample.put(nonz_idx[sample_idx],0)
                else:
                    # Delete Time
                    nonz_idx = time.nonzero()[0]
                    time.put(nonz_idx[time_idx],0)
            
            self._compute_MSR(time, gene, sample)

        return time, gene, sample

    def _multiple_node_deletion(self, time, gene, sample):
        self._compute_MSR(time,gene,sample)

        if self.n_gene < 50 or self.n_sample <50 or self.n_time >= 50:
            pass
        else:  
            self._compute_MSR(time, gene, sample)

            while (self.MSR > self.delta):
                            
                deleted = 0

                gene_to_del = self.MSR_gene > (self.l * self.MSR)
                nonz_idx = gene.nonzero()[0]
                if (gene_to_del.any()):
                    deleted = 1
                gene.put(nonz_idx[gene_to_del], 0)
                self._compute_MSR(time, gene, sample)

                sample_to_del = self.MSR_sample > (self.l * self.MSR)
                nonz_idx = sample.nonzero()[0]
                if (sample_to_del.any()):
                    deleted = 1
                sample.put(nonz_idx[sample_to_del], 0)
                self._compute_MSR(time, gene, sample)

                time_to_del = self.MSR_time > (self.l * self.MSR)
                nonz_idx = samples.nonzero()[0]
                if (time_to_del.any()):
                    deleted = 1
                time.put(nonz_idx[time_to_del], 0)
                self._compute_MSR(time, gene, sample)

                if (not deleted):
                    break

 #           self._compute_MSR(time, gene, sample)

        return time, gene, sample

    def _node_addition(self, time, gene, sample):

        if self.i >1:
            self.D[self.t, self.g, self.s] = self.D_asli


        while True:
            
            self._compute_MSR(time, gene, sample)
            n_time = np.count_nonzero(time)
            n_gene = np.count_nonzero(gene)
            n_sample = np.count_nonzero(sample)

            # gene addition
            elems_to_add = self.MSR_gene <= self.MSR
            nonz_idx = gene.nonzero()[0]
            gene.put(nonz_idx[elems_to_add], 1)
            #recalculate MSR
            if elems_to_add.any():
                self._compute_MSR(time, gene, sample, gene_add=True)

            # sample addition
            elems_to_add = self.MSR_sample<= self.MSR
            nonz_idx = sample.nonzero()[0]
            sample.put(nonz_idx[elems_to_add], 1)
            #recalculate
            if elems_to_add.any():
                self._compute_MSR(time, gene, sample, sample_add=True)

            # time addition
            elems_to_add = self.MSR_time <= self.MSR
            nonz_idx = time.nonzero()[0]
            time.put(nonz_idx[elems_to_add], 1)

            if (n_gene== np.count_nonzero(gene)) and \
               (n_sample == np.count_nonzero(sample)) and \
               (n_time == np.count_nonzero(time)):
                break

        return time, gene, sample

    def _mask(self, time, gene, sample, minval, maxval):

        self.g = np.expand_dims(np.expand_dims(np.nonzero(gene)[0], axis=1), axis=0)
        self.s = np.expand_dims(np.expand_dims(np.nonzero(sample)[0], axis=0), axis=0)
        self.t = np.expand_dims(np.expand_dims(np.nonzero(time)[0], axis=1), axis=2)
        if (self.mask_mode == 'random'):
            #menyimpan nilai asli untuk dipanggil lagi pada node addition
            self.D_asli = D[self.t, self.g, self.s]
            # membuat nilai random pada D
            shape = np.count_nonzero(time), np.count_nonzero(gene), np.count_nonzero(sample)
            mask_vals = np.random.uniform(minval, maxval, shape)
            self.D[self.t, self.g, self.s] = mask_vals
        else:
            self.D[self.t, self.g, self.s] = np.nan

    def fit(self, delta=2.5, l=1.005, time_cutoff=50, gene_cutoff=50,
            sample_cutoff=50, tol=1e-5, mask_mode='nan', verbose=False):

        self.delta = delta
        self.l = l
        self.time_cutoff = time_cutoff
        self.gene_cutoff = gene_cutoff
        self.sample_cutoff = sample_cutoff
        self.tol = tol
        self.mask_mode = mask_mode

        n_time, n_gene, n_sample = self.D.shape
        minval, maxval = np.nanmin(self.D), np.nanmax(self.D)

        result_gene = []
        result_sample = []
        result_time = []

        self.i = 1       
        
        while True:
            
            if (verbose):
                print(self.i)
            time = np.ones(n_time, dtype=np.bool)
            gene = np.ones(n_gene, dtype=np.bool)
            sample = np.ones(n_sample, dtype=np.bool)

            # Multiple node deletion
            time, gene, sample = self._multiple_node_deletion(time,
                                                                  gene,
                                                                  sample)
           

            # Single node deletion
            time, gene, sample = self._single_node_deletion(time,
                                                                gene,
                                                                sample)
           

            # Node addition
            time, gene, sample = self._node_addition(time,
                                                         gene,
                                                         sample)            
            

            # Check for trivial tricluster
            if (time.sum() == 1) or (gene.sum() == 1) or (sample.sum() == 1):
                break  # trivial bicluster
            # Check if the aren't any unused values in D
            if ((mask_mode == 'nan') and (np.isnan(self.D).all())):
                break

            # Mask values
            self._mask(time, gene, sample, minval, maxval)

            result_time.append(time)
            result_gene.append(gene)
            result_sample.append(sample)
            if (verbose):
                print("--- MSR = " + str(self.MSR))

            self.i += 1

        self.result_time = result_time
        self.result_gene = result_gene
        self.result_sample = result_sample

    def get_triclusters(self):
        """
        Returns the triclusters found by the algorithm.
        """
        return self.result_time, self.result_gene, self.result_sample



