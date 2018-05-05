class Data_Directories:
    intervals = {
            'day0': '/mnt/lab_data/kundaje/projects/skin/data/bds/processed.atac.2016-10-09.bioreps_peaks_signals/primary_keratinocyte-d00.GGR.Stanford_Greenleaf.ATAC-seq/peak/idr/optimal_set/primary_keratinocyte-d00.GGR.Stanford_Greenleaf.ATAC-seq_rep1-rep2.IDR0.1.filt.narrowPeak.gz',
            'day3': '/mnt/lab_data/kundaje/projects/skin/data/bds/processed.atac.2016-10-09.bioreps_peaks_signals/primary_keratinocyte-d30.GGR.Stanford_Greenleaf.ATAC-seq/peak/idr/optimal_set/primary_keratinocyte-d30.GGR.Stanford_Greenleaf.ATAC-seq_ppr.IDR0.1.filt.narrowPeak.gz',
            'day6': '/mnt/lab_data/kundaje/projects/skin/data/bds/processed.atac.2016-10-09.bioreps_peaks_signals/primary_keratinocyte-d60.GGR.Stanford_Greenleaf.ATAC-seq/peak/idr/conservative_set/primary_keratinocyte-d60.GGR.Stanford_Greenleaf.ATAC-seq_rep1-rep2.IDR0.1.filt.narrowPeak.gz'
            }
    input_atac = {
            'day0': {
                '100': '/srv/scratch/csfoo/projects/ggr/bcolz/atac/Day-0.0-merged.trim.PE2SE.fixed.nodup.namesorted.fragments.bed.gz-100bp-even',
                '140': '/srv/scratch/csfoo/projects/ggr/bcolz/atac/Day-0.0-merged.trim.PE2SE.fixed.nodup.namesorted.fragments.bed.gz-140bp-nuc'
                },
            'day3': {
                '100': '/srv/scratch/csfoo/projects/ggr/bcolz/atac/Day-3.0-merged.trim.PE2SE.fixed.nodup.namesorted.fragments.bed.gz-100bp-even',
                '140': '/srv/scratch/csfoo/projects/ggr/bcolz/atac/Day-3.0-merged.trim.PE2SE.fixed.nodup.namesorted.fragments.bed.gz-140bp-nuc'
                },
            'day6':{
                '100': '/srv/scratch/csfoo/projects/ggr/bcolz/atac/Day-6.0-merged.trim.PE2SE.fixed.nodup.namesorted.fragments.bed.gz-100bp-even',
                '140': '/srv/scratch/csfoo/projects/ggr/bcolz/atac/Day-6.0-merged.trim.PE2SE.fixed.nodup.namesorted.fragments.bed.gz-140bp-nuc'
                }
            }
    output_histone = {
            'day0':{
                'H3K27ac': '/mnt/lab_data/kundaje/projects/skin/data/bds/processed.chipseq.2017-01-23.histones/primary_keratinocyte-d0.GGR.Stanford_Snyder.ChIP-seq_H3K27ac/signal/macs2/pooled_rep/primary_keratinocyte-d0.GGR.Stanford_Snyder.ChIP-seq_H3K27ac.b1.t1.ACAGTG.1.PE2SE.nodup_pooled.tagAlign_x_ChIP-seq_input_H3K27ac.b1.1.PE2SE.nodup_pooled.tagAlign.pval.signal.bw',
                'H3K27me3': '/mnt/lab_data/kundaje/projects/skin/data/bds/processed.chipseq.2017-01-23.histones/primary_keratinocyte-d0.GGR.Stanford_Snyder.ChIP-seq_H3K27me3/signal/macs2/pooled_rep/primary_keratinocyte-d0.GGR.Stanford_Snyder.ChIP-seq_H3K27me3.b1.t1.ATCACG.1.PE2SE.nodup_pooled.tagAlign_x_ChIP-seq_input_H3K27me3.b1.1.PE2SE.nodup_pooled.tagAlign.pval.signal.bw',
                'H3K4me1': '/mnt/lab_data/kundaje/projects/skin/data/bds/processed.chipseq.2017-01-23.histones/primary_keratinocyte-d0.GGR.Stanford_Snyder.ChIP-seq_H3K4me1/signal/macs2/pooled_rep/primary_keratinocyte-d0.GGR.Stanford_Snyder.ChIP-seq_H3K4me1.b1.t1.ACAGTG.1.PE2SE.nodup_pooled.tagAlign_x_ChIP-seq_input_H3K4me1.b1.1.PE2SE.nodup_pooled.tagAlign.pval.signal.bw'
                },
            'day3':{
                'H3K27ac': '/mnt/lab_data/kundaje/projects/skin/data/bds/processed.chipseq.2017-01-23.histones/primary_keratinocyte-d3.GGR.Stanford_Snyder.ChIP-seq_H3K27ac/signal/macs2/pooled_rep/primary_keratinocyte-d3.GGR.Stanford_Snyder.ChIP-seq_H3K27ac.b1.t1.ACTTGA.1.PE2SE.nodup_pooled.tagAlign_x_ChIP-seq_input_H3K27ac.b1.1.PE2SE.nodup_pooled.tagAlign.pval.signal.bw',
                'H3K27me3': '/mnt/lab_data/kundaje/projects/skin/data/bds/processed.chipseq.2017-01-23.histones/primary_keratinocyte-d3.GGR.Stanford_Snyder.ChIP-seq_H3K27me3/signal/macs2/pooled_rep/primary_keratinocyte-d3.GGR.Stanford_Snyder.ChIP-seq_H3K27me3.b1.t1.TTAGGC.1.PE2SE.nodup_pooled.tagAlign_x_ChIP-seq_input_H3K27me3.b1.1.PE2SE.nodup_pooled.tagAlign.pval.signal.bw',
                'H3K4me1': '/mnt/lab_data/kundaje/projects/skin/data/bds/processed.chipseq.2017-01-23.histones/primary_keratinocyte-d3.GGR.Stanford_Snyder.ChIP-seq_H3K4me1/signal/macs2/pooled_rep/primary_keratinocyte-d3.GGR.Stanford_Snyder.ChIP-seq_H3K4me1.b1.t1.ACTTGA.1.PE2SE.nodup_pooled.tagAlign_x_ChIP-seq_input_H3K4me1.b1.1.PE2SE.nodup_pooled.tagAlign.pval.signal.bw'
                },
            'day6':{
                'H3K27ac': '/mnt/lab_data/kundaje/projects/skin/data/bds/processed.chipseq.2017-01-23.histones/primary_keratinocyte-d6.GGR.Stanford_Snyder.ChIP-seq_H3K27ac/signal/macs2/pooled_rep/primary_keratinocyte-d6.GGR.Stanford_Snyder.ChIP-seq_H3K27ac.b1.t1.GGCTAC.1.PE2SE.nodup_pooled.tagAlign_x_ChIP-seq_input_H3K27ac.b1.1.PE2SE.nodup_pooled.tagAlign.pval.signal.bw',
                'H3K27me3': '/mnt/lab_data/kundaje/projects/skin/data/bds/processed.chipseq.2017-01-23.histones/primary_keratinocyte-d6.GGR.Stanford_Snyder.ChIP-seq_H3K27me3/signal/macs2/pooled_rep/primary_keratinocyte-d6.GGR.Stanford_Snyder.ChIP-seq_H3K27me3.b1.t1.CAGATC.1.PE2SE.nodup_pooled.tagAlign_x_ChIP-seq_input_H3K27me3.b1.1.PE2SE.nodup_pooled.tagAlign.pval.signal.bw',
                'H3K4me1': '/mnt/lab_data/kundaje/projects/skin/data/bds/processed.chipseq.2017-01-23.histones/primary_keratinocyte-d6.GGR.Stanford_Snyder.ChIP-seq_H3K4me1/signal/macs2/pooled_rep/primary_keratinocyte-d6.GGR.Stanford_Snyder.ChIP-seq_H3K4me1.b1.t1.GGCTAC.1.PE2SE.nodup_pooled.tagAlign_x_ChIP-seq_input_H3K4me1.b1.1.PE2SE.nodup_pooled.tagAlign.pval.signal.bw'
                }
            }
