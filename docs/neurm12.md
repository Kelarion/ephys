---
title: cerebellum project
---

# Supplementary material

These are the files and figures which supplement the full report on the deep cerebellar nuclei. See [this github folder] for the code to generate the figures from the report. (See 'main_analysis.m' -- it may require the rest of the ephys repository to run.)

[this github folder]:https://github.com/Kelarion/ephys/tree/master/cerebellum

Movies
------
**Sup. Mov. 1**: [muad_dib_6-19-17_2217_meanEventOnset.avi]

> This dataset had the clearest oscillatory events, and thus we are more confident about the detected onset times. Only cells with significant phase-locked spiking were considered for this analysis. 

![alt text][screencap]

[muad_dib_6-19-17_2217_meanEventOnset.avi]: https://drive.google.com/file/d/1ZMPyG3y3KX1GGi99VKt-0TK4YMBdshu4/view?usp=sharing

[screencap]:img/meanMovieStill.jpg "A yellow dot calms the thousand-limbed beast"

Online experimental methods
------
> These are the experimental steps I am familiar with, which could have an impact on the results of the report.

**Electrophysiology.**

Following one or two days habituation to the ladder, craniotomies were made over the cerebellar cortex to allow access to the deep nuclei. These did occasionally damage the tissue, which might contribute to the phenomena described in the report. In the first recordings (Muad'Dib), a NeuroNexus 4x16 probe was used for acquisition; in subsequent recordings (all other mice), a custom design (Hires 4x16 Janelia probes) was used. In a recording session, we began acquiring  recordings once online phototagging confirmed DCN localization (see *Probe localization*); after ~10 minutes, probes were lowered at least 300 microns until the next acquisition to prevent overlap in neurons sampled in each recording. In some recordings, recordings of spontaneous activity (i.e. without overt movement) were also taken -- however, these are yet unprocessed and thus were not analysed in the report.



**Spike sorting.**

Software developed on-site (JRClust, James Jun) was used to isolate individual neurons from the extracellular recording. A manual curation step was used to assess automatic performance and label well-isolated units from multiunit activity. I was conservative in my assessment, and only well-isolated units were included in the analysis.



**Probe localization.**

Blue light was shined on the surface of the cerebellum to excited ChR2-expressin Purkinje cells. Light was delivered in volleys of short pulses (each around 2ms), with a variable number of pulses in each volley. The duration of pulses and frequency of their presentation were varied in some recordings, but again these are not yet spike-sorted and could not be included. In addition to online photostimulation, probes were localized post-mortem by visualization of DiI tracts *in vitro*.

