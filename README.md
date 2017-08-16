A piece of code for fitting infrared continuum measurements to spectral energy distribution templates. So far, I have implemented the Kirkpatrick et al. ([2012](http://adsabs.harvard.edu/abs/2012ApJ...759..139K); [2015](http://adsabs.harvard.edu/abs/2015ApJ...814....9K)) SED libraries. The code accepts measured fluxes from *Herschel* PACS and ALMA.

In the `data` directory I have included Kirkpatrick et al. (2012; 2015) templates (private communication) as well as filter throughputs from the *Herschel* instruments (which can be found [here](http://svo2.cab.inta-csic.es/svo/theory/fps/index.php?mode=browse&gname=Herschel)). For reference and my future convenience, here are the [PACS documentation](http://svo2.cab.inta-csic.es/svo/theory/fps/index.php?mode=browse&gname=Herschel) and [SPIRE
handbook](http://herschel.esac.esa.int/Docs/SPIRE/html/spire_om.html).

Here's what the data tree looks like:

    data
    ├── herschel_filters
    │   ├── Herschel-Pacs.blue.dat
    │   ├── Herschel-Pacs.green.dat
    │   ├── Herschel-Pacs.red.dat
    │   ├── Herschel-SPIRE.PLW.dat
    │   ├── Herschel-SPIRE.PLW_ext.dat
    │   ├── Herschel-SPIRE.PMW.dat
    │   ├── Herschel-SPIRE.PMW_ext.dat
    │   ├── Herschel-SPIRE.PSW.dat
    │   └── Herschel-SPIRE.PSW_ext.dat
    ├── kirkpatrick+12
    │   ├── featureless_AGN_SED.txt
    │   ├── silicate_AGN_SED.txt
    │   ├── z1_SF_SED.txt
    │   └── z2_SF_SED.txt
    └── kirkpatrick+15
        ├── Color_based_library
        │   ├── IR_COLOR1.txt
        │   ├── IR_COLOR2.txt
        │   ├── IR_COLOR3.txt
        │   ├── IR_COLOR4.txt
        │   ├── IR_COLOR5.txt
        │   ├── IR_COLOR6.txt
        │   ├── IR_COLOR7.txt
        │   └── IR_COLOR8.txt
        ├── Comprehensive_library
        │   ├── AGN1.txt
        │   ├── AGN2.txt
        │   ├── AGN3.txt
        │   ├── AGN4.txt
        │   ├── Composite1.txt
        │   ├── Composite2.txt
        │   ├── Composite3.txt
        │   ├── Composite4.txt
        │   ├── SFG1.txt
        │   ├── SFG2.txt
        │   └── SFG3.txt
        └── MIR_based_library
            ├── MIR0.0.txt
            ├── MIR0.1.txt
            ├── MIR0.2.txt
            ├── MIR0.3.txt
            ├── MIR0.4.txt
            ├── MIR0.5.txt
            ├── MIR0.6.txt
            ├── MIR0.7.txt
            ├── MIR0.8.txt
            ├── MIR0.9.txt
            └── MIR1.0.txt
