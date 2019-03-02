/*	VISIT: Vegetation Integrative SImulator for Trace gases				*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & RIGC/JAMSTEC						*/
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version 1.0.0	cerated in August 14, 2007							*/

/* History */
/* Created 2007 / 12 /25				*/
/* Revised 2008 / 02 /17				*/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include"structure.h"
#include"prototype.h"

/* open input files **************************************************/
void open_input(
	FILE *fp_s[IFILEN], 
	FILE *fp_c[4]
){
    char filename[128], telumoid[16];

    if(ISIMIP_RUN==0){
        if(DL_HCLIM == 102){
            /* UEA-CRU data from 1901 - 2002 */
            if( (fp_c[0]=fopen("./data/cru21_cld_1901-2002.dat","rt"))==NULL ){  
                printf("No cru21_cld_1901-2002.dat\n");  
                exit(1); 
            }
            if( (fp_c[1]=fopen("./data/cru21_pre_1901-2002.dat","rt"))==NULL ){  
                printf("No cru21_pre_1901-2002.dat\n");  
                exit(1); 
            }
            if( (fp_c[2]=fopen("./data/cru21_tmp_1901-2002.dat","rt"))==NULL ){  
                printf("No cru21_tmp_1901-2002.dat\n");  
                exit(1); 
            }
            if( (fp_c[3]=fopen("./data/cru21_vap_1901-2002.dat","rt"))==NULL ){  
                printf("No cru21_vap_1901-2002.dat\n");  
                exit(1); 
            }
        }else if(DL_HCLIM==106){
            /* UEA-CRU data from 1901 - 2006: 2010/01/04 (A.Ito) */
            if( (fp_c[0]=fopen("./data/cru_ts30_1901-2006.cld.dat","rt"))==NULL ){  
                printf("No cru_ts30_1901-2006.cld.dat\n");  
                exit(1); 
            }
            if( (fp_c[1]=fopen("./data/cru_ts30_1901-2006.pre.dat","rt"))==NULL ){  
                printf("No cru_ts30_1901-2006.pre.dat\n");  
                exit(1); 
            }
            if( (fp_c[2]=fopen("./data/cru_ts30_1901-2006.tmp.dat","rt"))==NULL ){  
                printf("No cru_ts30_1901-2006.tmp.dat\n");  
                exit(1); 
            }
            if( (fp_c[3]=fopen("./data/cru_ts30_1901-2006.vap.dat","rt"))==NULL ){  
                printf("No cru_ts30_1901-2006.vap.dat\n");  
                exit(1); 
            }
        }else if(DL_HCLIM==109){
            /* UEA-CRU data from 1901 - 2009: 2011/05/02 (A.Ito) */
            if( (fp_c[0]=fopen("./data/cru31_cld_1901-2009.txt","rt"))==NULL ){  
                printf("No cru31_cld_1901-2009.txt\n");  
                exit(1); 
            }
            if( (fp_c[1]=fopen("./data/cru31_pre_1901-2009.txt","rt"))==NULL ){  
                printf("No cru31_pre_1901-2009.txt\n");  
                exit(1); 
            }
            if( (fp_c[2]=fopen("./data/cru31_tmp_1901-2009.txt","rt"))==NULL ){  
                printf("No cru31_tmp_1901-2009.txt\n");  
                exit(1); 
            }
            if( (fp_c[3]=fopen("./data/cru31_vap_1901-2009.txt","rt"))==NULL ){  
                printf("No cru31_vap_1901-2009.txt\n");  
                exit(1); 
            }
        }else if(DL_HCLIM==111){
            /* UEA-CRU data from 1901 - 2011: 2012/04/14 (A.Ito) */
            if( (fp_c[0]=fopen("./data/cru32_cld_1901-2011.txt","rt"))==NULL ){  
                printf("No cru32_cld_1901-2011.txt\n");  
                exit(1); 
            }
            if( (fp_c[1]=fopen("./data/cru32_pre_1901-2011.txt","rt"))==NULL ){  
                printf("No cru32_pre_1901-2011.txt\n");  
                exit(1); 
            }
            if( (fp_c[2]=fopen("./data/cru32_tmp_1901-2011.txt","rt"))==NULL ){  
                printf("No cru32_tmp_1901-2011.txt\n");  
                exit(1); 
            }
            if( (fp_c[3]=fopen("./data/cru32_vap_1901-2011.txt","rt"))==NULL ){  
                printf("No cru32_vap_1901-2011.txt\n");  
                exit(1); 
            }
        }else if(DL_HCLIM==112){
            /* UEA-CRU data from 1901 - 2012: 2013/11/26 (A.Ito) */
            if( (fp_c[0]=fopen("./data/cru321_cld_1901-2012.txt","rt"))==NULL ){  
                printf("No cru321_cld_1901-2012.txt\n");  
                exit(1); 
            }
            if( (fp_c[1]=fopen("./data/cru321_pre_1901-2012.txt","rt"))==NULL ){  
                printf("No cru321_pre_1901-2012.txt\n");  
                exit(1); 
            }
            if( (fp_c[2]=fopen("./data/cru321_tmp_1901-2012.txt","rt"))==NULL ){  
                printf("No cru321_tmp_1901-2012.txt\n");  
                exit(1); 
            }
            if( (fp_c[3]=fopen("./data/cru321_vap_1901-2012.txt","rt"))==NULL ){  
                printf("No cru321_vap_1901-2012.txt\n");  
                exit(1); 
            }
        }else if(DL_HCLIM==113){
            /* UEA-CRU data from 1901 - 2013: 2014/08/26 (A.Ito) */
            if( (fp_c[0]=fopen("./data/cru322_cld_1901-2013.txt","rt"))==NULL ){
                printf("No cru322_cld_1901-2013.txt\n");
                exit(1); 
            }
            if( (fp_c[1]=fopen("./data/cru322_pre_1901-2013.txt","rt"))==NULL ){
                printf("No cru322_pre_1901-2013.txt\n");
                exit(1); 
            }
            if( (fp_c[2]=fopen("./data/cru322_tmp_1901-2013.txt","rt"))==NULL ){
                printf("No cru322_tmp_1901-2013.txt\n");
                exit(1); 
            }
            if( (fp_c[3]=fopen("./data/cru322_vap_1901-2013.txt","rt"))==NULL ){
                printf("No cru322_vap_1901-2013.txt\n");
                exit(1); 
            }
        }else if(DL_HCLIM==114){
            /* UEA-CRU data from 1901 - 2014: 2015/11/10 (A.Ito) */
            if( (fp_c[0]=fopen("./data/cru323_cld_1901-2014.txt","rt"))==NULL ){
                printf("No cru323_cld_1901-2014.txt\n");
                exit(1); 
            }
            if( (fp_c[1]=fopen("./data/cru323_pre_1901-2014.txt","rt"))==NULL ){
                printf("No cru323_pre_1901-2014.txt\n");
                exit(1); 
            }
            if( (fp_c[2]=fopen("./data/cru323_tmp_1901-2014.txt","rt"))==NULL ){
                printf("No cru323_tmp_1901-2014.txt\n");
                exit(1); 
            }
            if( (fp_c[3]=fopen("./data/cru323_vap_1901-2014.txt","rt"))==NULL ){
                printf("No cru323_vap_1901-2014.txt\n");
                exit(1); 
            }
        }else if(DL_HCLIM==115){
            /* UEA-CRU data from 1901 - 2015: 2016/10/11 (A.Ito) */
            if( (fp_c[0]=fopen("./data/cru324_cld_1901-2015.txt","rt"))==NULL ){
                printf("No cru324_cld_1901-2014.txt\n");
                exit(1); 
            }
            if( (fp_c[1]=fopen("./data/cru324_pre_1901-2015.txt","rt"))==NULL ){
                printf("No cru324_pre_1901-2014.txt\n");
                exit(1); 
            }
            if( (fp_c[2]=fopen("./data/cru324_tmp_1901-2015.txt","rt"))==NULL ){
                printf("No cru324_tmp_1901-2014.txt\n");
                exit(1); 
            }
            if( (fp_c[3]=fopen("./data/cru324_vap_1901-2015.txt","rt"))==NULL ){
                printf("No cru324_vap_1901-2014.txt\n");
                exit(1); 
            }
        }else if(DL_HCLIM==116){
            /* UEA-CRU data from 1901 - 2016: 2016/10/18 (A.Ito) */
            if( (fp_c[0]=fopen("./data/cru325_cld_1901-2016.txt","rt"))==NULL ){
                printf("No cru325_cld_1901-2016.txt\n");
                exit(1);
            }
            if( (fp_c[1]=fopen("./data/cru325_pre_1901-2016.txt","rt"))==NULL ){
                printf("No cru325_pre_1901-2016.txt\n");
                exit(1);
            }
            if( (fp_c[2]=fopen("./data/cru325_tmp_1901-2016.txt","rt"))==NULL ){
                printf("No cru325_tmp_1901-2016.txt\n");
                exit(1);
            }
            if( (fp_c[3]=fopen("./data/cru325_vap_1901-2016.txt","rt"))==NULL ){
                printf("No cru325_vap_1901-2016.txt\n");
                exit(1);
            }
        }else if(DL_HCLIM==117){
            /* UEA-CRU data from 1901 - 2017: 2018/08/28 (A.Ito) */
            if( (fp_c[0]=fopen("./data/cru326_cld_1901-2017.txt","rt"))==NULL ){
                printf("No cru326_cld_1901-2017.txt\n");
                exit(1);
            }
            if( (fp_c[1]=fopen("./data/cru326_pre_1901-2017.txt","rt"))==NULL ){
                printf("No cru326_pre_1901-2017.txt\n");
                exit(1);
            }
            if( (fp_c[2]=fopen("./data/cru326_tmp_1901-2017.txt","rt"))==NULL ){
                printf("No cru326_tmp_1901-2017.txt\n");
                exit(1);
            }
            if( (fp_c[3]=fopen("./data/cru326_vap_1901-2017.txt","rt"))==NULL ){
                printf("No cru326_vap_1901-2017.txt\n");
                exit(1);
            }
        }else{
            printf("No CRU data\n");
            exit(1);
        }
    }else /* ISI-MIP: 2012/06/27 by A.Ito */
    /* ICARUS SSPs: 2016/08/14 by A.Ito */
    if(ISIMIP_RUN == 1){
        switch(SCENARIO_ID){
            case 2001: case 2005:
                case 2221: case 2223: case 2228: case 2230:
                case 2601: case 2602: case 2603:
                if( (fp_c[0]=fopen("./data/tas_rcp2p6_isimip_hadgem_1950-2099.flt","rb"))==NULL ){
                    printf("No tas_rcp2p6_isimip_hadgem_1950-2099.flt\n");  
                    exit(1); 
                }
                if( (fp_c[1]=fopen("./data/pr_rcp2p6_isimip_hadgem_1950-2099.flt","rb"))==NULL ){  
                    printf("No pr_rcp2p6_isimip_hadgem_1950-2099.flt\n");  
                    exit(1); 
                }
                if( (fp_c[2]=fopen("./data/rhs_rcp2p6_isimip_hadgem_1950-2099.flt","rb"))==NULL ){  
                    printf("No rhs_rcp2p6_isimip_hadgem_1950-2099.flt\n");  
                    exit(1); 
                }
                if( (fp_c[3]=fopen("./data/rsds_rcp2p6_isimip_hadgem_1950-2099.flt","rb"))==NULL ){  
                    printf("No rsds_rcp2p6_isimip_hadgem_1950-2099.flt\n");  
                    exit(1); 
                }
                break;
            case 2002: case 2006:
                if( (fp_c[0]=fopen("./data/tas_rcp8p5_isimip_hadgem_1950-2099.flt","rb"))==NULL ){  
                    printf("No tas_rcp8p5_isimip_hadgem_1950-2099.flt\n");  
                    exit(1); 
                }
                if( (fp_c[1]=fopen("./data/pr_rcp8p5_isimip_hadgem_1950-2099.flt","rb"))==NULL ){  
                    printf("No pr_rcp8p5_isimip_hadgem_1950-2099.flt\n");  
                    exit(1); 
                }
                if( (fp_c[2]=fopen("./data/rhs_rcp8p5_isimip_hadgem_1950-2099.flt","rb"))==NULL ){  
                    printf("No rhs_rcp8p5_isimip_hadgem_1950-2099.flt\n");  
                    exit(1); 
                }
                if( (fp_c[3]=fopen("./data/rsds_rcp8p5_isimip_hadgem_1950-2099.flt","rb"))==NULL ){  
                    printf("No rsds_rcp8p5_isimip_hadgem_1950-2099.flt\n");  
                    exit(1); 
                }
                break;
            case 2003: case 2007:
                case 2222: case 2224: case 2226: case 2229: case 2231:
                if( (fp_c[0]=fopen("./data/tas_rcp4p5_isimip_hadgem_1950-2099.flt","rb"))==NULL ){  
                    printf("No tas_rcp4p5_isimip_hadgem_1950-2099.flt\n");  
                    exit(1); 
                }
                if( (fp_c[1]=fopen("./data/pr_rcp4p5_isimip_hadgem_1950-2099.flt","rb"))==NULL ){  
                    printf("No pr_rcp4p5_isimip_hadgem_1950-2099.flt\n");  
                    exit(1); 
                }
                if( (fp_c[2]=fopen("./data/rhs_rcp4p5_isimip_hadgem_1950-2099.flt","rb"))==NULL ){  
                    printf("No rhs_rcp4p5_isimip_hadgem_1950-2099.flt\n");  
                    exit(1); 
                }
                if( (fp_c[3]=fopen("./data/rsds_rcp4p5_isimip_hadgem_1950-2099.flt","rb"))==NULL ){  
                    printf("No rsds_rcp4p5_isimip_hadgem_1950-2099.flt\n");  
                    exit(1); 
                }
                break;
            case 2004: case 2008:
                case 2225: case 2227:
                if( (fp_c[0]=fopen("./data/tas_rcp6p0_isimip_hadgem_1950-2099.flt","rb"))==NULL ){  
                    printf("No tas_rcp6p0_isimip_hadgem_1950-2099.flt\n");  
                    exit(1); 
                }
                if( (fp_c[1]=fopen("./data/pr_rcp6p0_isimip_hadgem_1950-2099.flt","rb"))==NULL ){  
                    printf("No pr_rcp6p0_isimip_hadgem_1950-2099.flt\n");  
                    exit(1); 
                }
                if( (fp_c[2]=fopen("./data/rhs_rcp6p0_isimip_hadgem_1950-2099.flt","rb"))==NULL ){  
                    printf("No rhs_rcp6p0_isimip_hadgem_1950-2099.flt\n");  
                    exit(1); 
                }
                if( (fp_c[3]=fopen("./data/rsds_rcp6p0_isimip_hadgem_1950-2099.flt","rb"))==NULL ){  
                    printf("No rsds_rcp6p0_isimip_hadgem_1950-2099.flt\n");  
                    exit(1); 
                }
                break;
            case 2011: case 2015:
                case 2241: case 2243: case 2248: case 2250:
                if( (fp_c[0]=fopen("./data/tas_rcp2p6_isimip_ipsl_1950-2099.flt","rb"))==NULL ){  
                    printf("No tas_rcp2p6_isimip_ipsl_1950-2099.flt\n");  
                    exit(1); 
                }
                if( (fp_c[1]=fopen("./data/pr_rcp2p6_isimip_ipsl_1950-2099.flt","rb"))==NULL ){  
                    printf("No pr_rcp2p6_isimip_ipsl_1950-2099.flt\n");  
                    exit(1); 
                }
                if( (fp_c[2]=fopen("./data/rhs_rcp2p6_isimip_ipsl_1950-2099.flt","rb"))==NULL ){  
                    printf("No huss_rcp2p6_isimip_ipsl_1950-2099.flt\n");  
                    exit(1); 
                }
                if( (fp_c[3]=fopen("./data/rsds_rcp2p6_isimip_ipsl_1950-2099.flt","rb"))==NULL ){  
                    printf("No rsds_rcp2p6_isimip_ipsl_1950-2099.flt\n");  
                    exit(1); 
                }
                break;
            case 2012: case 2016:
                if( (fp_c[0]=fopen("./data/tas_rcp8p5_isimip_ipsl_1950-2099.flt","rb"))==NULL ){  
                    printf("No tas_rcp8p5_isimip_ipsl_1950-2099.flt\n");  
                    exit(1); 
                }
                if( (fp_c[1]=fopen("./data/pr_rcp8p5_isimip_ipsl_1950-2099.flt","rb"))==NULL ){  
                    printf("No pr_rcp8p5_isimip_ipsl_1950-2099.flt\n");  
                    exit(1); 
                }
                if( (fp_c[2]=fopen("./data/rhs_rcp8p5_isimip_ipsl_1950-2099.flt","rb"))==NULL ){  
                    printf("No rhs_rcp8p5_isimip_ipsl_1950-2099.flt\n");  
                    exit(1); 
                }
                if( (fp_c[3]=fopen("./data/rsds_rcp8p5_isimip_ipsl_1950-2099.flt","rb"))==NULL ){  
                    printf("No rsds_rcp8p5_isimip_ipsl_1950-2099.flt\n");  
                    exit(1); 
                }
                break;
            case 2013: case 2017:
                case 2242: case 2244: case 2246: case 2249: case 2251:
                if( (fp_c[0]=fopen("./data/tas_rcp4p5_isimip_ipsl_1950-2099.flt","rb"))==NULL ){  
                    printf("No tas_rcp4p5_isimip_ipsl_1950-2099.flt\n");  
                    exit(1); 
                }
                if( (fp_c[1]=fopen("./data/pr_rcp4p5_isimip_ipsl_1950-2099.flt","rb"))==NULL ){  
                    printf("No pr_rcp4p5_isimip_ipsl_1950-2099.flt\n");  
                    exit(1); 
                }
                if( (fp_c[2]=fopen("./data/rhs_rcp4p5_isimip_ipsl_1950-2099.flt","rb"))==NULL ){  
                    printf("No rhs_rcp4p5_isimip_ipsl_1950-2099.flt\n");  
                    exit(1); 
                }
                if( (fp_c[3]=fopen("./data/rsds_rcp4p5_isimip_ipsl_1950-2099.flt","rb"))==NULL ){  
                    printf("No rsds_rcp4p5_isimip_ipsl_1950-2099.flt\n");  
                    exit(1); 
                }
                break;
            case 2014: case 2018:
                case 2245: case 2247:
                if( (fp_c[0]=fopen("./data/tas_rcp6p0_isimip_ipsl_1950-2099.flt","rb"))==NULL ){  
                    printf("No tas_rcp6p0_isimip_ipsl_1950-2099.flt\n");  
                    exit(1); 
                }
                if( (fp_c[1]=fopen("./data/pr_rcp6p0_isimip_ipsl_1950-2099.flt","rb"))==NULL ){  
                    printf("No pr_rcp6p0_isimip_ipsl_1950-2099.flt\n");  
                    exit(1); 
                }
                if( (fp_c[2]=fopen("./data/rhs_rcp6p0_isimip_ipsl_1950-2099.flt","rb"))==NULL ){  
                    printf("No rhs_rcp6p0_isimip_ipsl_1950-2099.flt\n");  
                    exit(1); 
                }
                if( (fp_c[3]=fopen("./data/rsds_rcp6p0_isimip_ipsl_1950-2099.flt","rb"))==NULL ){  
                    printf("No rsds_rcp6p0_isimip_ipsl_1950-2099.flt\n");  
                    exit(1); 
                }
                break;
            case 2021: case 2025:
                case 2201: case 2203: case 2208: case 2210:
                if( (fp_c[0]=fopen("./data/tas_rcp2p6_isimip_gfdl_1950-2099.flt","rb"))==NULL ){  
                    printf("No tas_rcp2p6_isimip_gfdl_1950-2099.flt\n");  
                    exit(1); 
                }
                if( (fp_c[1]=fopen("./data/pr_rcp2p6_isimip_gfdl_1950-2099.flt","rb"))==NULL ){  
                    printf("No pr_rcp2p6_isimip_gfdl_1950-2099.flt\n");  
                    exit(1); 
                }
                if( (fp_c[2]=fopen("./data/rhs_rcp2p6_isimip_gfdl_1950-2099.flt","rb"))==NULL ){  
                    printf("No rhs_rcp2p6_isimip_gfdl_1950-2099.flt\n");  
                    exit(1); 
                }
                if( (fp_c[3]=fopen("./data/rsds_rcp2p6_isimip_gfdl_1950-2099.flt","rb"))==NULL ){  
                    printf("No rsds_rcp2p6_isimip_gfdl_1950-2099.flt\n");  
                    exit(1); 
                }
                break;
            case 2022: case 2026:
                if( (fp_c[0]=fopen("./data/tas_rcp8p5_isimip_gfdl_1950-2099.flt","rb"))==NULL ){  
                    printf("No tas_rcp8p5_isimip_gfdl_1950-2099.flt\n");  
                    exit(1); 
                }
                if( (fp_c[1]=fopen("./data/pr_rcp8p5_isimip_gfdl_1950-2099.flt","rb"))==NULL ){  
                    printf("No pr_rcp8p5_isimip_gfdl_1950-2099.flt\n");  
                    exit(1); 
                }
                if( (fp_c[2]=fopen("./data/rhs_rcp8p5_isimip_gfdl_1950-2099.flt","rb"))==NULL ){  
                    printf("No rhs_rcp8p5_isimip_gfdl_1950-2099.flt\n");  
                    exit(1); 
                }
                if( (fp_c[3]=fopen("./data/rsds_rcp8p5_isimip_gfdl_1950-2099.flt","rb"))==NULL ){  
                    printf("No rsds_rcp8p5_isimip_gfdl_1950-2099.flt\n");  
                    exit(1); 
                }
                break;
            case 2023: case 2027:
                case 2202: case 2204: case 2206: case 2209: case 2211:
                if( (fp_c[0]=fopen("./data/tas_rcp4p5_isimip_gfdl_1950-2099.flt","rb"))==NULL ){  
                    printf("No tas_rcp4p5_isimip_gfdl_1950-2099.flt\n");  
                    exit(1); 
                }
                if( (fp_c[1]=fopen("./data/pr_rcp4p5_isimip_gfdl_1950-2099.flt","rb"))==NULL ){  
                    printf("No pr_rcp4p5_isimip_gfdl_1950-2099.flt\n");  
                    exit(1); 
                }
                if( (fp_c[2]=fopen("./data/rhs_rcp4p5_isimip_gfdl_1950-2099.flt","rb"))==NULL ){  
                    printf("No rhs_rcp4p5_isimip_gfdl_1950-2099.flt\n");  
                    exit(1); 
                }
                if( (fp_c[3]=fopen("./data/rsds_rcp4p5_isimip_gfdl_1950-2099.flt","rb"))==NULL ){  
                    printf("No rsds_rcp4p5_isimip_gfdl_1950-2099.flt\n");  
                    exit(1); 
                }
                break;
            case 2024: case 2028:
                case 2205: case 2207:
                if( (fp_c[0]=fopen("./data/tas_rcp6p0_isimip_gfdl_1950-2099.flt","rb"))==NULL ){  
                    printf("No tas_rcp6p0_isimip_gfdl_1950-2099.flt\n");  
                    exit(1); 
                }
                if( (fp_c[1]=fopen("./data/pr_rcp6p0_isimip_gfdl_1950-2099.flt","rb"))==NULL ){  
                    printf("No pr_rcp6p0_isimip_gfdl_1950-2099.flt\n");  
                    exit(1); 
                }
                if( (fp_c[2]=fopen("./data/rhs_rcp6p0_isimip_gfdl_1950-2099.flt","rb"))==NULL ){  
                    printf("No rhs_rcp6p0_isimip_gfdl_1950-2099.flt\n");  
                    exit(1); 
                }
                if( (fp_c[3]=fopen("./data/rsds_rcp6p0_isimip_gfdl_1950-2099.flt","rb"))==NULL ){  
                    printf("No rsds_rcp6p0_isimip_gfdl_1950-2099.flt\n");  
                    exit(1); 
                }
                break;
            case 2031: case 2035:
                case 2261: case 2263: case 2268: case 2270:
                if( (fp_c[0]=fopen("./data/tas_rcp2p6_isimip_miroc_1950-2099.flt","rb"))==NULL ){  
                    printf("No tas_rcp2p6_isimip_miroc_1950-2099.flt\n");  
                    exit(1); 
                }
                if( (fp_c[1]=fopen("./data/pr_rcp2p6_isimip_miroc_1950-2099.flt","rb"))==NULL ){  
                    printf("No pr_rcp2p6_isimip_miroc_1950-2099.flt\n");  
                    exit(1); 
                }
                if( (fp_c[2]=fopen("./data/rhs_rcp2p6_isimip_miroc_1950-2099.flt","rb"))==NULL ){  
                    printf("No rhs_rcp2p6_isimip_miroc_1950-2099.flt\n");  
                    exit(1); 
                }
                if( (fp_c[3]=fopen("./data/rsds_rcp2p6_isimip_miroc_1950-2099.flt","rb"))==NULL ){  
                    printf("No rsds_rcp2p6_isimip_miroc_1950-2099.flt\n");  
                    exit(1); 
                }
                break;
            case 2032: case 2036:
                if( (fp_c[0]=fopen("./data/tas_rcp8p5_isimip_miroc_1950-2099.flt","rb"))==NULL ){  
                    printf("No tas_rcp8p5_isimip_miroc_1950-2099.flt\n");  
                    exit(1); 
                }
                if( (fp_c[1]=fopen("./data/pr_rcp8p5_isimip_miroc_1950-2099.flt","rb"))==NULL ){  
                    printf("No pr_rcp8p5_isimip_miroc_1950-2099.flt\n");  
                    exit(1); 
                }
                if( (fp_c[2]=fopen("./data/rhs_rcp8p5_isimip_miroc_1950-2099.flt","rb"))==NULL ){  
                    printf("No rhs_rcp8p5_isimip_miroc_1950-2099.flt\n");  
                    exit(1); 
                }
                if( (fp_c[3]=fopen("./data/rsds_rcp8p5_isimip_miroc_1950-2099.flt","rb"))==NULL ){  
                    printf("No rsds_rcp8p5_isimip_miroc_1950-2099.flt\n");  
                    exit(1); 
                }
                break;
            case 2033: case 2037:
                case 2262: case 2264: case 2266: case 2269: case 2271:
                if( (fp_c[0]=fopen("./data/tas_rcp4p5_isimip_miroc_1950-2099.flt","rb"))==NULL ){  
                    printf("No tas_rcp4p5_isimip_miroc_1950-2099.flt\n");  
                    exit(1); 
                }
                if( (fp_c[1]=fopen("./data/pr_rcp4p5_isimip_miroc_1950-2099.flt","rb"))==NULL ){  
                    printf("No pr_rcp4p5_isimip_miroc_1950-2099.flt\n");  
                    exit(1); 
                }
                if( (fp_c[2]=fopen("./data/rhs_rcp4p5_isimip_miroc_1950-2099.flt","rb"))==NULL ){  
                    printf("No rhs_rcp4p5_isimip_miroc_1950-2099.flt\n");  
                    exit(1); 
                }
                if( (fp_c[3]=fopen("./data/rsds_rcp4p5_isimip_miroc_1950-2099.flt","rb"))==NULL ){  
                    printf("No rsds_rcp4p5_isimip_miroc_1950-2099.flt\n");  
                    exit(1); 
                }
                break;
            case 2034: case 2038:
                case 2265: case 2267:
                if( (fp_c[0]=fopen("./data/tas_rcp6p0_isimip_miroc_1950-2099.flt","rb"))==NULL ){  
                    printf("No tas_rcp6p0_isimip_miroc_1950-2099.flt\n");  
                    exit(1); 
                }
                if( (fp_c[1]=fopen("./data/pr_rcp6p0_isimip_miroc_1950-2099.flt","rb"))==NULL ){  
                    printf("No pr_rcp6p0_isimip_miroc_1950-2099.flt\n");  
                    exit(1); 
                }
                if( (fp_c[2]=fopen("./data/rhs_rcp6p0_isimip_miroc_1950-2099.flt","rb"))==NULL ){  
                    printf("No rhs_rcp6p0_isimip_miroc_1950-2099.flt\n");  
                    exit(1); 
                }
                if( (fp_c[3]=fopen("./data/rsds_rcp6p0_isimip_miroc_1950-2099.flt","rb"))==NULL ){  
                    printf("No rsds_rcp6p0_isimip_miroc_1950-2099.flt\n");  
                    exit(1); 
                }
                break;
            case 2041: case 2045:
                case 2281: case 2283: case 2288: case 2290:
                if( (fp_c[0]=fopen("./data/tas_rcp2p6_isimip_noresm_1950-2099.flt","rb"))==NULL ){  
                    printf("No tas_rcp2p6_isimip_noresm_1950-2099.flt\n");  
                    exit(1); 
                }
                if( (fp_c[1]=fopen("./data/pr_rcp2p6_isimip_noresm_1950-2099.flt","rb"))==NULL ){  
                    printf("No pr_rcp2p6_isimip_noresm_1950-2099.flt\n");  
                    exit(1); 
                }
                if( (fp_c[2]=fopen("./data/rhs_rcp2p6_isimip_noresm_1950-2099.flt","rb"))==NULL ){  
                    printf("No rhs_rcp2p6_isimip_noresm_1950-2099.flt\n");  
                    exit(1); 
                }
                if( (fp_c[3]=fopen("./data/rsds_rcp2p6_isimip_noresm_1950-2099.flt","rb"))==NULL ){  
                    printf("No rsds_rcp2p6_isimip_noresm_1950-2099.flt\n");  
                    exit(1); 
                }
                break;
            case 2042: case 2046:
                if( (fp_c[0]=fopen("./data/tas_rcp8p5_isimip_noresm_1950-2099.flt","rb"))==NULL ){  
                    printf("No tas_rcp8p5_isimip_noresm_1950-2099.flt\n");  
                    exit(1); 
                }
                if( (fp_c[1]=fopen("./data/pr_rcp8p5_isimip_noresm_1950-2099.flt","rb"))==NULL ){  
                    printf("No pr_rcp8p5_isimip_noresm_1950-2099.flt\n");  
                    exit(1); 
                }
                if( (fp_c[2]=fopen("./data/rhs_rcp8p5_isimip_noresm_1950-2099.flt","rb"))==NULL ){  
                    printf("No rhs_rcp8p5_isimip_noresm_1950-2099.flt\n");  
                    exit(1); 
                }
                if( (fp_c[3]=fopen("./data/rsds_rcp8p5_isimip_noresm_1950-2099.flt","rb"))==NULL ){  
                    printf("No rsds_rcp8p5_isimip_noresm_1950-2099.flt\n");  
                    exit(1); 
                }
                break;
            case 2043: case 2047:
                case 2282: case 2284: case 2286: case 2289: case 2291:
                if( (fp_c[0]=fopen("./data/tas_rcp4p5_isimip_noresm_1950-2099.flt","rb"))==NULL ){  
                    printf("No tas_rcp4p5_isimip_noresm_1950-2099.flt\n");  
                    exit(1); 
                }
                if( (fp_c[1]=fopen("./data/pr_rcp4p5_isimip_noresm_1950-2099.flt","rb"))==NULL ){  
                    printf("No pr_rcp4p5_isimip_noresm_1950-2099.flt\n");  
                    exit(1); 
                }
                if( (fp_c[2]=fopen("./data/rhs_rcp4p5_isimip_noresm_1950-2099.flt","rb"))==NULL ){  
                    printf("No rhs_rcp4p5_isimip_noresm_1950-2099.flt\n");  
                    exit(1); 
                }
                if( (fp_c[3]=fopen("./data/rsds_rcp4p5_isimip_noresm_1950-2099.flt","rb"))==NULL ){  
                    printf("No rsds_rcp4p5_isimip_noresm_1950-2099.flt\n");  
                    exit(1); 
                }
                break;
            case 2044: case 2048:
                case 2285: case 2287:
                if( (fp_c[0]=fopen("./data/tas_rcp6p0_isimip_noresm_1950-2099.flt","rb"))==NULL ){  
                    printf("No tas_rcp6p0_isimip_noresm_1950-2099.flt\n");  
                    exit(1); 
                }
                if( (fp_c[1]=fopen("./data/pr_rcp6p0_isimip_noresm_1950-2099.flt","rb"))==NULL ){  
                    printf("No pr_rcp6p0_isimip_noresm_1950-2099.flt\n");  
                    exit(1); 
                }
                if( (fp_c[2]=fopen("./data/rhs_rcp6p0_isimip_noresm_1950-2099.flt","rb"))==NULL ){  
                    printf("No rhs_rcp6p0_isimip_noresm_1950-2099.flt\n");  
                    exit(1); 
                }
                if( (fp_c[3]=fopen("./data/rsds_rcp6p0_isimip_noresm_1950-2099.flt","rb"))==NULL ){  
                    printf("No rsds_rcp6p0_isimip_noresm_1950-2099.flt\n");  
                    exit(1); 
                }
                break;
            default:
                printf("No ISI-MIP scenario\n");
                exit(1);
                break;
        }
    }else if(ISIMIP_RUN == 2){
        switch(SCENARIO_ID){
            case 4011: case 4012:
                if( (fp_c[0]=fopen("./data/tas_bced_1960_1999_gfdl-esm2m_spin_1901-1930_hist_1860-2005.flt","rb"))==NULL ){
                    printf("No tas_bced_1960_1999_gfdl-esm2m_spin_1901-1930_hist_1860-2005.flt\n");
                    exit(1); 
                }
                if( (fp_c[1]=fopen("./data/pr_bced_1960_1999_gfdl-esm2m_hist_spin_1901-1930_hist_1860-2005.flt","rb"))==NULL ){
                    printf("No pr_bced_1960_1999_gfdl-esm2m_hist_spin_1901-1930_hist_1860-2005.flt\n");
                    exit(1); 
                }
                if( (fp_c[2]=fopen("./data/hurs_bced_1960_1999_gfdl-esm2m_spin_1901-1930_hist_1901-2005.flt","rb"))==NULL ){
                    printf("No hurs_bced_1960_1999_gfdl-esm2m_spin_1901-1930_hist_1901-2005.flt\n");
                    exit(1); 
                }
                if( (fp_c[3]=fopen("./data/rsds_bced_1960_1999_gfdl-esm2m_spin_1901-1930_hist_1901-2005.flt","rb"))==NULL ){
                    printf("No rsds_bced_1960_1999_gfdl-esm2m_spin_1901-1930_hist_1901-2005.flt\n");
                    exit(1); 
                }
                break;
            case 4021: case 4022: case 4023: case 4024:
                if( (fp_c[0]=fopen("./data/tas_bced_1960_1999_ipsl-cm5a-lr_spin_1901-1930_hist_1901-2005.flt","rb"))==NULL ){
                    printf("No tas_bced_1960_1999_ipsl-cm5a-lr_spin_1901-1930_hist_1901-2005.flt\n");
                    exit(1); 
                }
                if( (fp_c[1]=fopen("./data/pr_bced_1960_1999_ipsl-cm5a-lr_spin_1901-1930_hist_1901-2005.flt","rb"))==NULL ){
                    printf("No pr_bced_1960_1999_ipsl-cm5a-lr_spin_1901-1930_hist_1901-2005.flt\n");
                    exit(1); 
                }
                if( (fp_c[2]=fopen("./data/hurs_bced_1960_1999_ipsl-cm5a-lr_spin_1901-1930_hist_1901-2005.flt","rb"))==NULL ){
                    printf("No hurs_bced_1960_1999_ipsl-cm5a-lr_spin_1901-1930_hist_1901-2005.flt\n");
                    exit(1); 
                }
                if( (fp_c[3]=fopen("./data/rsds_bced_1960_1999_ipsl-cm5a-lr_spin_1901-1930_hist_1901-2005.flt","rb"))==NULL ){
                    printf("No rsds_bced_1960_1999_ipsl-cm5a-lr_spin_1901-1930_hist_1901-2005.flt\n");
                    exit(1); 
                }
                break;
        default:
            printf("No ISI-MIP scenario\n");
            exit(1);
            break;
        }
    }else if(ISIMIP_RUN == 3){
        switch(SCENARIO_ID){
            case 5001:
                if( (fp_c[0]=fopen("./data/tas_gswp3_1901_1930_1901_2010_ver2.flt","rb"))==NULL ){
                    printf("No tas_gswp3_1901_1930_1901_2010_ver2.flt\n");
                    exit(1); 
                }
                if( (fp_c[1]=fopen("./data/pr_gswp3_1901_1930_1901_2010_ver2.flt","rb"))==NULL ){
                    printf("No pr_gswp3_1901_1930_1901_2010_ver2.flt\n");
                    exit(1); 
                }
                if( (fp_c[2]=fopen("./data/huss_gswp3_1901_1930_1901_2010_ver2.flt","rb"))==NULL ){
                    printf("No huss_gswp3_1901_1930_1901_2010_ver2.flt\n");
                    exit(1); 
                }
                if( (fp_c[3]=fopen("./data/rsds_gswp3_1901_1930_1901_2010_ver2.flt","rb"))==NULL ){
                    printf("No rsds_gswp3_1901_1930_1901_2010_ver2.flt\n");
                    exit(1); 
                }
                break;
             case 5002: /* corrected (pgfv2.1: 2010 => 2012): 2016/05/06 by A.Ito */
                if( (fp_c[0]=fopen("./data/tas_pgfv2.1_1901_1930_1901_2012_ver2.flt","rb"))==NULL ){
                    printf("No tas_pgfv2.1_1901_1930_1901_2012_ver2.flt\n");
                    exit(1); 
                }
                if( (fp_c[1]=fopen("./data/pr_pgfv2.1_1901_1930_1901_2012_ver2.flt","rb"))==NULL ){
                    printf("No pr_pgfv2.1_1901_1930_1901_2012_ver2.flt\n");
                    exit(1); 
                }
                if( (fp_c[2]=fopen("./data/huss_pgfv2.1_1901_1930_1901_2012_ver2.flt","rb"))==NULL ){
                    printf("No huss_pgfv2.1_1901_1930_1901_2012_ver2.flt\n");
                    exit(1); 
                }
                if( (fp_c[3]=fopen("./data/rsds_pgfv2.1_1901_1930_1901_2012_ver2.flt","rb"))==NULL ){
                    printf("No rsds_pgfv2.1_1901_1930_1901_2012_ver2.flt\n");
                    exit(1); 
                }
                break;
             case 5003:
                if( (fp_c[0]=fopen("./data/tas_watch_1901_1930_1901_2001_ver2.flt","rb"))==NULL ){
                    printf("No tas_watch_1901_1930_1901_2001_ver2.flt\n");
                    exit(1); 
                }
                if( (fp_c[1]=fopen("./data/pr_watch_1901_1930_1901_2001_ver2.flt","rb"))==NULL ){
                    printf("No pr_watch_1901_1930_1901_2001_ver2.flt\n");
                    exit(1); 
                }
                if( (fp_c[2]=fopen("./data/huss_watch_1901_1930_1901_2001_ver2.flt","rb"))==NULL ){
                    printf("No huss_watch_1901_1930_1901_2001_ver2.flt\n");
                    exit(1); 
                }
                if( (fp_c[3]=fopen("./data/rsds_watch_1901_1930_1901_2001_ver2.flt","rb"))==NULL ){
                    printf("No rsds_watch_1901_1930_1901_2001_ver2.flt\n");
                    exit(1); 
                }
                break;
             case 5004:
                if( (fp_c[0]=fopen("./data/tas_watch+wfdei_1901_1930_1901_2010_ver2.flt","rb"))==NULL ){
                    printf("No tas_watch+wfdei_1901_1930_1901_2010_ver2.flt\n");
                    exit(1); 
                }
                if( (fp_c[1]=fopen("./data/pr_watch+wfdei_1901_1930_1901_2010_ver2.flt","rb"))==NULL ){
                    printf("No pr_watch+wfdei_1901_1930_1901_2010_ver2.flt\n");
                    exit(1); 
                }
                if( (fp_c[2]=fopen("./data/huss_watch+wfdei_1901_1930_1901_2010_ver2.flt","rb"))==NULL ){
                    printf("No huss_watch+wfdei_1901_1930_1901_2010_ver2.flt\n");
                    exit(1); 
                }
                if( (fp_c[3]=fopen("./data/rsds_watch+wfdei_1901_1930_1901_2010_ver2.flt","rb"))==NULL ){
                    printf("No rsds_watch+wfdei_1901_1930_1901_2010_ver2.flt\n");
                    exit(1); 
                }
                break;
        }
    }else if(ISIMIP_RUN == 4){
        /* ISI-MIP2.1b: 2016/12/27 by A.Ito */
        switch(SCENARIO_ID){
            case 5010: case 5013: case 5014: case 5060:
                if( (fp_c[0]=fopen("./data/tas_mon_GFDL-ESM2M_r1i1p1_EWEMBI_landonly_pi_1661-2299.flt","rb"))==NULL ){
                    printf("No tas_mon_GFDL-ESM2M_r1i1p1_EWEMBI_landonly_pi_1661-2299.flt\n");
                    exit(1); 
                }
                if( (fp_c[1]=fopen("./data/pr_mon_GFDL-ESM2M_r1i1p1_EWEMBI_landonly_pi_1661-2299.flt","rb"))==NULL ){
                    printf("No pr_mon_GFDL-ESM2M_r1i1p1_EWEMBI_landonly_pi_1661-2299.flt\n");
                    exit(1); 
                }
                if( (fp_c[2]=fopen("./data/huss_mon_GFDL-ESM2M_r1i1p1_EWEMBI_landonly_pi_1661-2299.flt","rb"))==NULL ){
                    printf("No huss_mon_GFDL-ESM2M_r1i1p1_EWEMBI_landonly_pi_1661-2299.flt\n");
                    exit(1); 
                }
                if( (fp_c[3]=fopen("./data/rsds_mon_GFDL-ESM2M_r1i1p1_EWEMBI_landonly_pi_1661-2299.flt","rb"))==NULL ){
                    printf("No rsds_mon_GFDL-ESM2M_r1i1p1_EWEMBI_landonly_pi_1661-2299.flt\n");
                    exit(1); 
                }
                break;
            case 5011: case 5015: case 5061:
                if( (fp_c[0]=fopen("./data/tas_mon_GFDL-ESM2M_r1i1p1_EWEMBI_landonly_pihir26_1661-2299.flt","rb"))==NULL ){
                    printf("No tas_mon_GFDL-ESM2M_r1i1p1_EWEMBI_landonly_pihir26_1661-2299.flt\n");
                    exit(1); 
                }
                if( (fp_c[1]=fopen("./data/pr_mon_GFDL-ESM2M_r1i1p1_EWEMBI_landonly_pihir26_1661-2299.flt","rb"))==NULL ){
                    printf("No pr_mon_GFDL-ESM2M_r1i1p1_EWEMBI_landonly_pihir26_1661-2299.flt\n");
                    exit(1); 
                }
                if( (fp_c[2]=fopen("./data/huss_mon_GFDL-ESM2M_r1i1p1_EWEMBI_landonly_pihir26_1661-2299.flt","rb"))==NULL ){
                    printf("No huss_mon_GFDL-ESM2M_r1i1p1_EWEMBI_landonly_pihir26_1661-2299.flt\n");
                    exit(1); 
                }
                if( (fp_c[3]=fopen("./data/rsds_mon_GFDL-ESM2M_r1i1p1_EWEMBI_landonly_pihir26_1661-2299.flt","rb"))==NULL ){
                    printf("No rsds_mon_GFDL-ESM2M_r1i1p1_EWEMBI_landonly_pihir26_1661-2299.flt\n");
                    exit(1); 
                }
                break;
            case 5012: case 5016: case 5062: case 5063:
                if( (fp_c[0]=fopen("./data/tas_mon_GFDL-ESM2M_r1i1p1_EWEMBI_landonly_pihir60_1661-2299.flt","rb"))==NULL ){
                    printf("No tas_mon_GFDL-ESM2M_r1i1p1_EWEMBI_landonly_pi_1661-2299.flt\n");
                    exit(1); 
                }
                if( (fp_c[1]=fopen("./data/pr_mon_GFDL-ESM2M_r1i1p1_EWEMBI_landonly_pihir60_1661-2299.flt","rb"))==NULL ){
                    printf("No pr_mon_GFDL-ESM2M_r1i1p1_EWEMBI_landonly_pihir60_1661-2299.flt\n");
                    exit(1); 
                }
                if( (fp_c[2]=fopen("./data/huss_mon_GFDL-ESM2M_r1i1p1_EWEMBI_landonly_pihir60_1661-2299.flt","rb"))==NULL ){
                    printf("No huss_mon_GFDL-ESM2M_r1i1p1_EWEMBI_landonly_pihir60_1661-2299.flt\n");
                    exit(1); 
                }
                if( (fp_c[3]=fopen("./data/rsds_mon_GFDL-ESM2M_r1i1p1_EWEMBI_landonly_pihir60_1661-2299.flt","rb"))==NULL ){
                    printf("No rsds_mon_GFDL-ESM2M_r1i1p1_EWEMBI_landonly_pihir60_1661-2299.flt\n");
                    exit(1); 
                }
                break;
            case 5020: case 5023: case 5024: case 5070:
                if( (fp_c[0]=fopen("./data/tas_mon_IPSL-CM5A-LR_r1i1p1_EWEMBI_landonly_pi_1661-2299.flt","rb"))==NULL ){
                    printf("No tas_mon_IPSL-CM5A-LR_r1i1p1_EWEMBI_landonly_pi_1661-2299.flt\n");
                    exit(1); 
                }
                if( (fp_c[1]=fopen("./data/pr_mon_IPSL-CM5A-LR_r1i1p1_EWEMBI_landonly_pi_1661-2299.flt","rb"))==NULL ){
                    printf("No pr_mon_IPSL-CM5A-LR_r1i1p1_EWEMBI_landonly_pi_1661-2299.flt\n");
                    exit(1); 
                }
                if( (fp_c[2]=fopen("./data/huss_mon_IPSL-CM5A-LR_r1i1p1_EWEMBI_landonly_pi_1661-2299.flt","rb"))==NULL ){
                    printf("No huss_mon_IPSL-CM5A-LR_r1i1p1_EWEMBI_landonly_pi_1661-2299.flt\n");
                    exit(1); 
                }
                if( (fp_c[3]=fopen("./data/rsds_mon_IPSL-CM5A-LR_r1i1p1_EWEMBI_landonly_pi_1661-2299.flt","rb"))==NULL ){
                    printf("No rsds_mon_IPSL-CM5A-LR_r1i1p1_EWEMBI_landonly_pi_1661-2299.flt\n");
                    exit(1); 
                }
                break;
            case 5021: case 5025: case 5071:
                if( (fp_c[0]=fopen("./data/tas_mon_IPSL-CM5A-LR_r1i1p1_EWEMBI_landonly_pihir26_1661-2299.flt","rb"))==NULL ){
                    printf("No tas_mon_IPSL-CM5A-LR_r1i1p1_EWEMBI_landonly_pihir26_1661-2299.flt\n");
                    exit(1); 
                }
                if( (fp_c[1]=fopen("./data/pr_mon_IPSL-CM5A-LR_r1i1p1_EWEMBI_landonly_pihir26_1661-2299.flt","rb"))==NULL ){
                    printf("No pr_mon_IPSL-CM5A-LR_r1i1p1_EWEMBI_landonly_pihir26_1661-2299.flt\n");
                    exit(1); 
                }
                if( (fp_c[2]=fopen("./data/huss_mon_IPSL-CM5A-LR_r1i1p1_EWEMBI_landonly_pihir26_1661-2299.flt","rb"))==NULL ){
                    printf("No huss_mon_IPSL-CM5A-LR_r1i1p1_EWEMBI_landonly_pihir26_1661-2299.flt\n");
                    exit(1); 
                }
                if( (fp_c[3]=fopen("./data/rsds_mon_IPSL-CM5A-LR_r1i1p1_EWEMBI_landonly_pihir26_1661-2299.flt","rb"))==NULL ){
                    printf("No rsds_mon_IPSL-CM5A-LR_r1i1p1_EWEMBI_landonly_pihir26_1661-2299.flt\n");
                    exit(1); 
                }
                break;
            case 5022: case 5026: case 5072: case 5073:
                if( (fp_c[0]=fopen("./data/tas_mon_IPSL-CM5A-LR_r1i1p1_EWEMBI_landonly_pihir60_1661-2299.flt","rb"))==NULL ){
                    printf("No tas_mon_IPSL-CM5A-LR_r1i1p1_EWEMBI_landonly_pihir60_1661-2299.flt\n");
                    exit(1); 
                }
                if( (fp_c[1]=fopen("./data/pr_mon_IPSL-CM5A-LR_r1i1p1_EWEMBI_landonly_pihir60_1661-2299.flt","rb"))==NULL ){
                    printf("No pr_mon_IPSL-CM5A-LR_r1i1p1_EWEMBI_landonly_pihir60_1661-2299.flt\n");
                    exit(1); 
                }
                if( (fp_c[2]=fopen("./data/huss_mon_IPSL-CM5A-LR_r1i1p1_EWEMBI_landonly_pihir60_1661-2299.flt","rb"))==NULL ){
                    printf("No huss_mon_IPSL-CM5A-LR_r1i1p1_EWEMBI_landonly_pihir60_1661-2299.flt\n");
                    exit(1); 
                }
                if( (fp_c[3]=fopen("./data/rsds_mon_IPSL-CM5A-LR_r1i1p1_EWEMBI_landonly_pihir60_1661-2299.flt","rb"))==NULL ){
                    printf("No rsds_mon_IPSL-CM5A-LR_r1i1p1_EWEMBI_landonly_pihir60_1661-2299.flt\n");
                    exit(1); 
                }
                break;
            case 5030: case 5033: case 5034: case 5080:
                if( (fp_c[0]=fopen("./data/tas_mon_MIROC5_r1i1p1_EWEMBI_landonly_pi_1661-2299.flt","rb"))==NULL ){
                    printf("No tas_mon_MIROC5_r1i1p1_EWEMBI_landonly_pi_1661-2299.flt\n");
                    exit(1); 
                }
                if( (fp_c[1]=fopen("./data/pr_mon_MIROC5_r1i1p1_EWEMBI_landonly_pi_1661-2299.flt","rb"))==NULL ){
                    printf("No pr_mon_MIROC5_r1i1p1_EWEMBI_landonly_pi_1661-2299.flt\n");
                    exit(1); 
                }
                if( (fp_c[2]=fopen("./data/huss_mon_MIROC5_r1i1p1_EWEMBI_landonly_pi_1661-2299.flt","rb"))==NULL ){
                    printf("No huss_mon_MIROC5_r1i1p1_EWEMBI_landonly_pi_1661-2299.flt\n");
                    exit(1); 
                }
                if( (fp_c[3]=fopen("./data/rsds_mon_MIROC5_r1i1p1_EWEMBI_landonly_pi_1661-2299.flt","rb"))==NULL ){
                    printf("No rsds_mon_MIROC5_r1i1p1_EWEMBI_landonly_pi_1661-2299.flt\n");
                    exit(1); 
                }
                break;
            case 5031: case 5035: case 5081:
                if( (fp_c[0]=fopen("./data/tas_mon_MIROC5_r1i1p1_EWEMBI_landonly_pihir26_1661-2299.flt","rb"))==NULL ){
                    printf("No tas_mon_MIROC5_r1i1p1_EWEMBI_landonly_pihir26_1661-2299.flt\n");
                    exit(1); 
                }
                if( (fp_c[1]=fopen("./data/pr_mon_MIROC5_r1i1p1_EWEMBI_landonly_pihir26_1661-2299.flt","rb"))==NULL ){
                    printf("No pr_mon_MIROC5_r1i1p1_EWEMBI_landonly_pihir26_1661-2299.flt\n");
                    exit(1); 
                }
                if( (fp_c[2]=fopen("./data/huss_mon_MIROC5_r1i1p1_EWEMBI_landonly_pihir26_1661-2299.flt","rb"))==NULL ){
                    printf("No huss_mon_MIROC5_r1i1p1_EWEMBI_landonly_pihir26_1661-2299.flt\n");
                    exit(1); 
                }
                if( (fp_c[3]=fopen("./data/rsds_mon_MIROC5_r1i1p1_EWEMBI_landonly_pihir26_1661-2299.flt","rb"))==NULL ){
                    printf("No rsds_mon_MIROC5_r1i1p1_EWEMBI_landonly_pihir26_1661-2299.flt\n");
                    exit(1); 
                }
                break;
            case 5032: case 5036: case 5082: case 5083:
                if( (fp_c[0]=fopen("./data/tas_mon_MIROC5_r1i1p1_EWEMBI_landonly_pihir60_1661-2299.flt","rb"))==NULL ){
                    printf("No tas_mon_MIROC5_r1i1p1_EWEMBI_landonly_pihir60_1661-2299.flt\n");
                    exit(1); 
                }
                if( (fp_c[1]=fopen("./data/pr_mon_MIROC5_r1i1p1_EWEMBI_landonly_pihir60_1661-2299.flt","rb"))==NULL ){
                    printf("No pr_mon_MIROC5_r1i1p1_EWEMBI_landonly_pihir60_1661-2299.flt\n");
                    exit(1); 
                }
                if( (fp_c[2]=fopen("./data/huss_mon_MIROC5_r1i1p1_EWEMBI_landonly_pihir60_1661-2299.flt","rb"))==NULL ){
                    printf("No huss_mon_MIROC5_r1i1p1_EWEMBI_landonly_pihir60_1661-2299.flt\n");
                    exit(1); 
                }
                if( (fp_c[3]=fopen("./data/rsds_mon_MIROC5_r1i1p1_EWEMBI_landonly_pihir60_1661-2299.flt","rb"))==NULL ){
                    printf("No rsds_mon_MIROC5_r1i1p1_EWEMBI_landonly_pihir60_1661-2299.flt\n");
                    exit(1); 
                }
                break;
            case 5040: case 5043: case 5044: case 5090:
                if( (fp_c[0]=fopen("./data/tas_mon_HadGEM2-ES_r1i1p1_EWEMBI_landonly_pi_1661-2299.flt","rb"))==NULL ){
                    printf("No tas_mon_HadGEM2-ES_r1i1p1_EWEMBI_landonly_pi_1661-2299.flt\n");
                    exit(1); 
                }
                if( (fp_c[1]=fopen("./data/pr_mon_HadGEM2-ES_r1i1p1_EWEMBI_landonly_pi_1661-2299.flt","rb"))==NULL ){
                    printf("No pr_mon_HadGEM2-ES_r1i1p1_EWEMBI_landonly_pi_1661-2299.flt\n");
                    exit(1); 
                }
                if( (fp_c[2]=fopen("./data/huss_mon_HadGEM2-ES_r1i1p1_EWEMBI_landonly_pi_1661-2299.flt","rb"))==NULL ){
                    printf("No huss_mon_HadGEM2-ES_r1i1p1_EWEMBI_landonly_pi_1661-2299.flt\n");
                    exit(1); 
                }
                if( (fp_c[3]=fopen("./data/rsds_mon_HadGEM2-ES_r1i1p1_EWEMBI_landonly_pi_1661-2299.flt","rb"))==NULL ){
                    printf("No rsds_mon_HadGEM2-ES_r1i1p1_EWEMBI_landonly_pi_1661-2299.flt\n");
                    exit(1); 
                }
                break;
            case 5041: case 5045: case 5091:
                if( (fp_c[0]=fopen("./data/tas_mon_HadGEM2-ES_r1i1p1_EWEMBI_landonly_pihir26_1661-2299.flt","rb"))==NULL ){
                    printf("No tas_mon_HadGEM2-ES_r1i1p1_EWEMBI_landonly_pihir26_1661-2299.flt\n");
                    exit(1); 
                }
                if( (fp_c[1]=fopen("./data/pr_mon_HadGEM2-ES_r1i1p1_EWEMBI_landonly_pihir26_1661-2299.flt","rb"))==NULL ){
                    printf("No pr_mon_HadGEM2-ES_r1i1p1_EWEMBI_landonly_pihir26_1661-2299.flt\n");
                    exit(1); 
                }
                if( (fp_c[2]=fopen("./data/huss_mon_HadGEM2-ES_r1i1p1_EWEMBI_landonly_pihir26_1661-2299.flt","rb"))==NULL ){
                    printf("No huss_mon_HadGEM2-ES_r1i1p1_EWEMBI_landonly_pihir26_1661-2299.flt\n");
                    exit(1); 
                }
                if( (fp_c[3]=fopen("./data/rsds_mon_HadGEM2-ES_r1i1p1_EWEMBI_landonly_pihir26_1661-2299.flt","rb"))==NULL ){
                    printf("No rsds_mon_HadGEM2-ES_r1i1p1_EWEMBI_landonly_pihir26_1661-2299.flt\n");
                    exit(1); 
                }
                break;
            case 5042: case 5046: case 5092: case 5093:
                if( (fp_c[0]=fopen("./data/tas_mon_HadGEM2-ES_r1i1p1_EWEMBI_landonly_pihir60_1661-2299.flt","rb"))==NULL ){
                    printf("No tas_mon_HadGEM2-ES_r1i1p1_EWEMBI_landonly_pihir60_1661-2299.flt\n");
                    exit(1); 
                }
                if( (fp_c[1]=fopen("./data/pr_mon_HadGEM2-ES_r1i1p1_EWEMBI_landonly_pihir60_1661-2299.flt","rb"))==NULL ){
                    printf("No pr_mon_HadGEM2-ES_r1i1p1_EWEMBI_landonly_pihir60_1661-2299.flt\n");
                    exit(1); 
                }
                if( (fp_c[2]=fopen("./data/huss_mon_HadGEM2-ES_r1i1p1_EWEMBI_landonly_pihir60_1661-2299.flt","rb"))==NULL ){
                    printf("No huss_mon_HadGEM2-ES_r1i1p1_EWEMBI_landonly_pihir60_1661-2299.flt\n");
                    exit(1); 
                }
                if( (fp_c[3]=fopen("./data/rsds_mon_HadGEM2-ES_r1i1p1_EWEMBI_landonly_pihir60_1661-2299.flt","rb"))==NULL ){
                    printf("No rsds_mon_HadGEM2-ES_r1i1p1_EWEMBI_landonly_pihir60_1661-2299.flt\n");
                    exit(1); 
                }
                break;
             case 5017:
                if( (fp_c[0]=fopen("./data/tas_mon_GFDL-ESM2M_r1i1p1_EWEMBI_landonly_pihir60_1661-2299.flt","rb"))==NULL ){
                    printf("No tas_mon_GFDL-ESM2M_r1i1p1_EWEMBI_landonly_pihir60_1661-2299.flt\n");
                    exit(1);
                }
                if( (fp_c[1]=fopen("./data/pr_mon_GFDL-ESM2M_r1i1p1_EWEMBI_landonly_pihir60_1661-2299.flt","rb"))==NULL ){
                    printf("No pr_mon_GFDL-ESM2M_r1i1p1_EWEMBI_landonly_pihir60_1661-2299.flt\n");
                    exit(1);
                }
                if( (fp_c[2]=fopen("./data/huss_mon_GFDL-ESM2M_r1i1p1_EWEMBI_landonly_pihir60_1661-2299.flt","rb"))==NULL ){
                    printf("No huss_mon_GFDL-ESM2M_r1i1p1_EWEMBI_landonly_pihir60_1661-2299.flt\n");
                    exit(1);
                }
                if( (fp_c[3]=fopen("./data/rsds_mon_GFDL-ESM2M_r1i1p1_EWEMBI_landonly_pihir60_1661-2299.flt","rb"))==NULL ){
                    printf("No rsds_mon_GFDL-ESM2M_r1i1p1_EWEMBI_landonly_pihir60_1661-2299.flt\n");
                    exit(1);
                }
                break;
             case 5027:
                if( (fp_c[0]=fopen("./data/tas_mon_IPSL-CM5A-LR_r1i1p1_EWEMBI_landonly_pihir60_1661-2299.flt","rb"))==NULL ){
                    printf("No tas_mon_IPSL-CM5A-LR_r1i1p1_EWEMBI_landonly_pihir60_1661-2299.flt\n");
                    exit(1);
                }
                if( (fp_c[1]=fopen("./data/pr_mon_IPSL-CM5A-LR_r1i1p1_EWEMBI_landonly_pihir60_1661-2299.flt","rb"))==NULL ){
                    printf("No pr_mon_IPSL-CM5A-LR_r1i1p1_EWEMBI_landonly_pihir60_1661-2299.flt\n");
                    exit(1);
                }
                if( (fp_c[2]=fopen("./data/huss_mon_IPSL-CM5A-LR_r1i1p1_EWEMBI_landonly_pihir60_1661-2299.flt","rb"))==NULL ){
                    printf("No huss_mon_IPSL-CM5A-LR_r1i1p1_EWEMBI_landonly_pihir60_1661-2299.flt\n");
                    exit(1);
                }
                if( (fp_c[3]=fopen("./data/rsds_mon_IPSL-CM5A-LR_r1i1p1_EWEMBI_landonly_pihir60_1661-2299.flt","rb"))==NULL ){
                    printf("No rsds_mon_IPSL-CM5A-LR_r1i1p1_EWEMBI_landonly_pihir60_1661-2299.flt\n");
                    exit(1);
                }
                break;
             case 5037:
                if( (fp_c[0]=fopen("./data/tas_mon_MIROC5_r1i1p1_EWEMBI_landonly_pihir60_1661-2299.flt","rb"))==NULL ){
                    printf("No tas_mon_MIROC5_r1i1p1_EWEMBI_landonly_pihir60_1661-2299.flt\n");
                    exit(1);
                }
                if( (fp_c[1]=fopen("./data/pr_mon_MIROC5_r1i1p1_EWEMBI_landonly_pihir60_1661-2299.flt","rb"))==NULL ){
                    printf("No pr_mon_MIROC5_r1i1p1_EWEMBI_landonly_pihir60_1661-2299.flt\n");
                    exit(1);
                }
                if( (fp_c[2]=fopen("./data/huss_mon_MIROC5_r1i1p1_EWEMBI_landonly_pihir60_1661-2299.flt","rb"))==NULL ){
                    printf("No huss_mon_MIROC5_r1i1p1_EWEMBI_landonly_pihir60_1661-2299.flt\n");
                    exit(1);
                }
                if( (fp_c[3]=fopen("./data/rsds_mon_MIROC5_r1i1p1_EWEMBI_landonly_pihir60_1661-2299.flt","rb"))==NULL ){
                    printf("No rsds_mon_MIROC5_r1i1p1_EWEMBI_landonly_pihir60_1661-2299.flt\n");
                    exit(1);
                }
                break;
             case 5047:
                if( (fp_c[0]=fopen("./data/tas_mon_HadGEM2-ES_r1i1p1_EWEMBI_landonly_pihir60_1661-2299.flt","rb"))==NULL ){
                    printf("No tas_mon_HadGEM2-ES_r1i1p1_EWEMBI_landonly_pihir60_1661-2299.flt\n");
                    exit(1);
                }
                if( (fp_c[1]=fopen("./data/pr_mon_HadGEM2-ES_r1i1p1_EWEMBI_landonly_pihir60_1661-2299.flt","rb"))==NULL ){
                    printf("No pr_mon_HadGEM2-ES_r1i1p1_EWEMBI_landonly_pihir60_1661-2299.flt\n");
                    exit(1);
                }
                if( (fp_c[2]=fopen("./data/huss_mon_HadGEM2-ES_r1i1p1_EWEMBI_landonly_pihir60_1661-2299.flt","rb"))==NULL ){
                    printf("No huss_mon_HadGEM2-ES_r1i1p1_EWEMBI_landonly_pihir60_1661-2299.flt\n");
                    exit(1);
                }
                if( (fp_c[3]=fopen("./data/rsds_mon_HadGEM2-ES_r1i1p1_EWEMBI_landonly_pihir60_1661-2299.flt","rb"))==NULL ){
                    printf("No rsds_mon_HadGEM2-ES_r1i1p1_EWEMBI_landonly_pihir60_1661-2299.flt\n");
                    exit(1);
                }
                break;
       }
    }
    
    /*********************************************************************/
	if( (fp_s[0]=fopen("./data/nat_reg_05.dat","rt"))==NULL ){  
		printf("No nat_reg_05.dat\n");  
		exit(1); 
	}else{
        Flag_FOPEN[0] ++;
    }
	if( (fp_s[1]=fopen("./data/olson_biome3_05.dat","rt"))==NULL ){
		printf("No olson_biome3_05.dat\n");  
		exit(1); 
	}else{
        Flag_FOPEN[1] ++;
    }
	if( (fp_s[2]=fopen("./data/climate_05.dat","rt"))==NULL ){
		printf("No climate_05.dat\n");  
		exit(1); 
	}else{
        Flag_FOPEN[2] ++;
    }
	if( (fp_s[3]=fopen("./data/geo_prop_05.dat","rt"))==NULL ){
		printf("No geo_prop_05.dat\n");  
		exit(1); 
	}else{
        Flag_FOPEN[3] ++;
    }
	if( (fp_s[4]=fopen("./data/prec_cru.dat","rt"))==NULL ){
		printf("No prec_cru.dat\n");  
		exit(1); 
	}else{
        Flag_FOPEN[4] ++;
    }
	if( (fp_s[5]=fopen("./data/rivnum_trip_05.dat","rt"))==NULL ){
		printf("No rivnum_trip_05.dat\n");  
		exit(1); 
	}else{
        Flag_FOPEN[5] ++;
    }
	if( (fp_s[6]=fopen("./data/Albedo_soil_EX.dat","rt"))==NULL ){
		printf("No Albedo_soil_EX.dat\n");  
		exit(1); 
	}else{
        Flag_FOPEN[6] ++;
    }
	if( (fp_s[7]=fopen("./data/pcnt_orgmatter.dat","rt"))==NULL ){  
		printf("No pcnt_orgmatter.dat\n");  
		exit(1); 
	}else{
        Flag_FOPEN[7] ++;
    }
	if( (fp_s[8]=fopen("./data/fls_slope.dat","rt"))==NULL ){  
		printf("No fls_slope.dat\n");  
		exit(1); 
	}else{
        Flag_FOPEN[8] ++;
    }
	if( (fp_s[9]=fopen("./data/fk_erodibility.dat","rt"))==NULL ){  
		printf("No fk_erodibility.dat\n");  
		exit(1); 
	}else{
        Flag_FOPEN[9] ++;
    }
	if( (fp_s[10]=fopen("./data/potveg_sage_05.dat","rt"))==NULL ){  
		printf("No potveg_sage_05.dat\n");  
		exit(1); 
	}else{
        Flag_FOPEN[10] ++;
    }
	if( (fp_s[11]=fopen("./data/luc_1700-1992.dat","rt"))==NULL ){  
		printf("No luc_1700-1992.dat\n");  
		exit(1); 
	}else{
        Flag_FOPEN[11] ++;
    }
	if( (fp_s[12]=fopen("./data/FAO_soil.textur.1nnegl.dat","rt"))==NULL ){  
		printf("No FAO_soil.textur.1nnegl.dat\n");  exit(1); 
	}else{
        Flag_FOPEN[12] ++;
    }
	if( (fp_s[13]=fopen("./data/bulkdens.dat","rt"))==NULL ){  
		printf("No bulkdens.dat\n");  
		exit(1); 
	}else{
        Flag_FOPEN[13] ++;
    }
	if( (fp_s[14]=fopen("./data/fieldcap.dat","rt"))==NULL ){  
		printf("No fieldcap.dat\n");  
		exit(1); 
	}else{
        Flag_FOPEN[14] ++;
    }
	if( (fp_s[15]=fopen("./data/soil_field_cap0-150_30min.dat","rt"))==NULL ){  
		printf("No soil_field_cap0-150_30min.dat\n");  
		exit(1); 
	}else{
        Flag_FOPEN[15] ++;
    }
	if( (fp_s[16]=fopen("./data/soil_field_cap0-30_30min.dat","rt"))==NULL ){  
		printf("No soil_field_cap0-30_30min.dat\n");  
		exit(1); 
	}else{
        Flag_FOPEN[16] ++;
    }
	if( (fp_s[17]=fopen("./data/soil_bulk_dens0-30_30min.dat","rt"))==NULL ){  
		printf("No soil_bulk_dens0-30_30min.dat\n");  
		exit(1); 
	}else{
        Flag_FOPEN[17] ++;
    }
	if( (fp_s[18]=fopen("./data/soil_clay_perc0-30_30min.dat","rt"))==NULL ){  
		printf("No soil_clay_perc0-30_30min.dat\n");  
		exit(1); 
	}else{
        Flag_FOPEN[18] ++;
    }
	if( (fp_s[19]=fopen("./data/soil_sand_perc0-30_30min.dat","rt"))==NULL ){  
		printf("No soil_sand_perc0-30_30min.dat\n");  
		exit(1); 
	}else{
        Flag_FOPEN[19] ++;
    }
	if( (fp_s[20]=fopen("./data/rastert_soilph1.dat","rt"))==NULL ){  
		printf("No rastert_soilph1.dat\n");  
		exit(1); 
	}else{
        Flag_FOPEN[20] ++;
    }

	if(ALT_FWETLAND == 1){
		/* Alternative data (NASA/GISS): 2011/03/30 by A.Ito */
		if( (fp_s[21]=fopen("./data/fwetland_giss.txt","rt"))==NULL ){  
			printf("No fwetland_giss.txt\n");  
			exit(1); 
		}else{
            Flag_FOPEN[21] ++;
        }
	}else{
		/* if( (fp_s[21]=fopen("./data/wetland_glwd_05.dat","rt"))==NULL ){  
		 printf("No wetland_glwd_05.dat\n");  
		 exit(1); 
		 } */
		/* revised wetland data: by A.Ito (2009/07/14) */
		if( (fp_s[21]=fopen("./data/wetland_glwd_31.dat","rt"))==NULL ){  
			printf("No wetland_glwd_31.dat\n");  
			exit(1); 
		}else{
            Flag_FOPEN[21] ++;
        }
	}
    /* WSL-Perogon welnad map: 2014/2/4 by A.Ito */
    if( (fp_s[83]=fopen("./data/fwetland_peregon.txt","rt"))==NULL ){
        printf("No fwetland_peregon.txt\n");  
        exit(1); 
    }else{
        Flag_FOPEN[83] ++;
    }
	
    /* GCP-CH4 inundation data: revised 2014/08/17 by A.Ito */
    if(ALT_INUND == 3){
        if( (fp_s[84]=fopen("./data/fw_00_swe_v1-2_smoothed_1999-2013.txt","rt"))==NULL ){
            printf("No GCP-CH4 inundation data\n");
            exit(1);
        }else{
            Flag_FOPEN[84] ++;
        }
    }else if(ALT_INUND == 4){
        if( (fp_s[84]=fopen("./data/fw_00_swe_v1-2_smoothed_1999-2013_norice.txt","rt"))==NULL ){
            printf("No GCP-CH4 inundation data\n");
            exit(1);
        }else{
            Flag_FOPEN[84] ++;
        }
    }else if(ALT_INUND == 5){
        if( (fp_s[84]=fopen("./data/fw_00_swe_v1-2_smoothed_1999-2013_norice-nowater.txt","rt"))==NULL ){
            printf("No GCP-CH4 inundation data\n");
            exit(1);
        }else{
            Flag_FOPEN[84] ++;
        }
    }else if(ALT_INUND == 6){
        if( (fp_s[84]=fopen("./data/fw_swamps-glwd_2000-2012.txt","rt"))==NULL ){
            printf("No fw_swamps-glwd_2000-2012.txt data\n");
            exit(1);
        }else{
            Flag_FOPEN[84] ++;
        }
    }else if(ALT_INUND == 7 || ALT_INUND == 8){ /* GCP-CH4 v2: 2018/08/28 by A.Ito */
        if( (fp_s[84]=fopen("./data/fw_gcp-ch4_wetlands_2000-2017.txt","rt"))==NULL ){
            printf("No fw_gcp-ch4_wetlands_2000-2017.txt data\n");
            exit(1);
        }else{
            Flag_FOPEN[84] ++;
        }
    }else{
        if( (fp_s[84]=fopen("./data/fw_swamp-biascor_1999-2013.txt","rt"))==NULL ){
            printf("No GCP-CH4 inundation data\n");
            exit(1);
        }else{
            Flag_FOPEN[84] ++;
        }
	}

	if(ALT_FWETLAND == 1){
		/* Alternative data (SAGE): 2011/03/30 by A.Ito */
		if( (fp_s[22]=fopen("./data/fpaddy_iis_sage.txt","rt"))==NULL ){  
			printf("No major_crps_rice_0.50x0.50.asc\n");  
			exit(1); 
		}else{
            Flag_FOPEN[22] ++;
        }
	}else{
		if( (fp_s[22]=fopen("./data/rice_sage_05.dat","rt"))==NULL ){  
			printf("No rice_sage_05.dat\n");  
			exit(1); 
		}else{
            Flag_FOPEN[22] ++;
        }
	}
	
	if( (fp_s[24]=fopen("./data/soil_n_1m.dat","rt"))==NULL ){  
		printf("No soil_n_1m.dat\n");  
		exit(1); 
	}else{
        Flag_FOPEN[24] ++;
    }
	
    /* GlobAlbedo *******************************/
    if( (fp_s[58]=fopen("./data/GlobAlbedo_av.flt","rb"))==NULL ){
        printf("No GlobAlbedo_av.flt\n");
        exit(1);
    }else{
        Flag_FOPEN[58] ++;
    }
    
	/* land-use: historical *************************************/
    /* SAGEHYDE: Hurtt et al. (2006) */
	if(LANDUSE == 6){
		if( (fp_s[26]=fopen("./data/EOS_SAGEHYDE_1D_crop.dat","rt"))==NULL ){  
			printf("No EOS_SAGEHYDE_1D_crop.dat\n");  
			exit(1); 
		}else{
            Flag_FOPEN[26] ++;
        }
		if( (fp_s[27]=fopen("./data/EOS_SAGEHYDE_1D_past.dat","rt"))==NULL ){  
			printf("No EOS_SAGEHYDE_1D_past.dat\n");  
			exit(1); 
		}else{
            Flag_FOPEN[27] ++;
        }
		if( (fp_s[28]=fopen("./data/EOS_SAGEHYDE_1D_prim.dat","rt"))==NULL ){  
			printf("No EOS_SAGEHYDE_1D_prim.dat\n");  
			exit(1); 
		}else{
            Flag_FOPEN[28] ++;
        }
		if( (fp_s[29]=fopen("./data/EOS_SAGEHYDE_1D_secd.dat","rt"))==NULL ){  
			printf("No EOS_SAGEHYDE_1D_secd.dat\n");  
			exit(1); 
		}else{
            Flag_FOPEN[29] ++;
        }
		if( (fp_s[30]=fopen("./data/EOS_SAGEHYDE_1D_ssma.dat","rt"))==NULL ){  
			printf("No EOS_SAGEHYDE_1D_ssma.dat\n");  
			exit(1); 
		}else{
            Flag_FOPEN[30] ++;
        }
		if( (fp_s[31]=fopen("./data/EOS_SAGEHYDE_1D_ssmb.dat","rt"))==NULL ){  
			printf("No EOS_SAGEHYDE_1D_ssmb.dat\n");  
			exit(1); 
		}else{
            Flag_FOPEN[31] ++;
        }
		if( (fp_s[32]=fopen("./data/EOS_SAGEHYDE_1D_t_cp.dat","rt"))==NULL ){  
			printf("No EOS_SAGEHYDE_1D_t_cp.dat\n");  
			exit(1); 
		}else{
            Flag_FOPEN[32] ++;
        }
		if( (fp_s[33]=fopen("./data/EOS_SAGEHYDE_1D_t_cs.dat","rt"))==NULL ){  
			printf("No EOS_SAGEHYDE_1D_t_cs.dat\n");  
			exit(1); 
		}else{
            Flag_FOPEN[33] ++;
        }
		if( (fp_s[34]=fopen("./data/EOS_SAGEHYDE_1D_t_pc.dat","rt"))==NULL ){  
			printf("No EOS_SAGEHYDE_1D_t_pc.dat\n");  
			exit(1); 
		}else{
            Flag_FOPEN[34] ++;
        }
		if( (fp_s[35]=fopen("./data/EOS_SAGEHYDE_1D_t_ps.dat","rt"))==NULL ){  
			printf("No EOS_SAGEHYDE_1D_t_ps.dat\n");  
			exit(1); 
		}else{
            Flag_FOPEN[35] ++;
        }
		if( (fp_s[36]=fopen("./data/EOS_SAGEHYDE_1D_t_sc.dat","rt"))==NULL ){  
			printf("No EOS_SAGEHYDE_1D_t_sc.dat\n");  
			exit(1); 
		}else{
            Flag_FOPEN[36] ++;
        }
		if( (fp_s[37]=fopen("./data/EOS_SAGEHYDE_1D_t_sp.dat","rt"))==NULL ){  
			printf("No EOS_SAGEHYDE_1D_t_sp.dat\n");  
			exit(1); 
		}else{
            Flag_FOPEN[37] ++;
        }
		if( (fp_s[38]=fopen("./data/EOS_SAGEHYDE_1D_t_ss1.dat","rt"))==NULL ){  
			printf("No EOS_SAGEHYDE_1D_t_ss1.dat\n");  
			exit(1); 
		}else{
            Flag_FOPEN[38] ++;
        }
		if( (fp_s[39]=fopen("./data/EOS_SAGEHYDE_1D_t_ss2.dat","rt"))==NULL ){  
			printf("No EOS_SAGEHYDE_1D_t_ss2.dat\n");  
			exit(1); 
		}else{
            Flag_FOPEN[39] ++;
        }
		if( (fp_s[40]=fopen("./data/EOS_SAGEHYDE_1D_t_ss3.dat","rt"))==NULL ){  
			printf("No EOS_SAGEHYDE_1D_t_ss3.dat\n");  
			exit(1); 
		}else{
            Flag_FOPEN[40] ++;
        }
		if( (fp_s[41]=fopen("./data/EOS_SAGEHYDE_1D_t_vc.dat","rt"))==NULL ){  
			printf("No EOS_SAGEHYDE_1D_t_vc.dat\n");  
			exit(1); 
		}else{
            Flag_FOPEN[41] ++;
        }
		if( (fp_s[42]=fopen("./data/EOS_SAGEHYDE_1D_t_vp.dat","rt"))==NULL ){  
			printf("No EOS_SAGEHYDE_1D_t_vp.dat\n");  
			exit(1); 
		}else{
            Flag_FOPEN[42] ++;
        }
		if( (fp_s[43]=fopen("./data/EOS_SAGEHYDE_1D_t_vs1.dat","rt"))==NULL ){  
			printf("No EOS_SAGEHYDE_1D_t_vs1.dat\n");  
			exit(1); 
		}else{
            Flag_FOPEN[43] ++;
        }
		if( (fp_s[44]=fopen("./data/EOS_SAGEHYDE_1D_t_vs2.dat","rt"))==NULL ){  
			printf("No EOS_SAGEHYDE_1D_t_vs2.dat\n");  
			exit(1); 
		}else{
            Flag_FOPEN[44] ++;
        }
	}
	/* U.NH harmonized historical land-use data (2010/01/31: by A.Ito) */
	/* URL  http://luh.unh.edu/ */
	if(LANDUSE == 8){
		if( (fp_s[26]=fopen("./data/luc_eos2_gcrop_1700-2005.dat","rt"))==NULL ){  
			printf("No luc_eos2_gcrop_1700-2005.dat\n");  
			exit(1); 
		}else{
            Flag_FOPEN[26] ++;
        }
		if( (fp_s[27]=fopen("./data/luc_eos2_gpast_1700-2005.dat","rt"))==NULL ){  
			printf("No luc_eos2_gpast_1700-2005.dat\n");  
			exit(1); 
		}else{
            Flag_FOPEN[27] ++;
        }
		if( (fp_s[28]=fopen("./data/luc_eos2_gothr_1700-2005.dat","rt"))==NULL ){  
			printf("No luc_eos2_gothr_1700-2005.dat\n");  
			exit(1); 
		}else{
            Flag_FOPEN[28] ++;
        }
		if( (fp_s[29]=fopen("./data/luc_eos2_gsecd_1700-2005.dat","rt"))==NULL ){  
			printf("No luc_eos2_gsecd_1700-2005.dat\n");  
			exit(1); 
		}else{
            Flag_FOPEN[29] ++;
        }
		if( (fp_s[30]=fopen("./data/luc_eos2_ssma_1700-2005.dat","rt"))==NULL ){  
			printf("No luc_eos2_ssma_1700-2005.dat\n");  
			exit(1); 
		}else{
            Flag_FOPEN[30] ++;
        }
		if( (fp_s[31]=fopen("./data/luc_eos2_ssmb_1700-2005.dat","rt"))==NULL ){  
			printf("No luc_eos2_ssmb_1700-2005.dat\n");  
			exit(1); 
		}else{
            Flag_FOPEN[31] ++;
        }
		if( (fp_s[32]=fopen("./data/luc_eos2_gflcp_1700-2005.dat","rt"))==NULL ){  
			printf("No luc_eos2_gflcp_1700-2005.dat\n");  
			exit(1); 
		}else{
            Flag_FOPEN[32] ++;
        }
		if( (fp_s[33]=fopen("./data/luc_eos2_gflcs_1700-2005.dat","rt"))==NULL ){  
			printf("No luc_eos2_gflcs_1700-2005.dat\n");  
			exit(1); 
		}else{
            Flag_FOPEN[33] ++;
        }
		if( (fp_s[34]=fopen("./data/luc_eos2_gflpc_1700-2005.dat","rt"))==NULL ){  
			printf("No luc_eos2_gflpc_1700-2005.dat\n");  
			exit(1); 
		}else{
            Flag_FOPEN[34] ++;
        }
		if( (fp_s[35]=fopen("./data/luc_eos2_gflps_1700-2005.dat","rt"))==NULL ){  
			printf("No luc_eos2_gflps_1700-2005.dat\n");  
			exit(1); 
		}else{
            Flag_FOPEN[35] ++;
        }
		if( (fp_s[36]=fopen("./data/luc_eos2_gflsc_1700-2005.dat","rt"))==NULL ){  
			printf("No luc_eos2_gflsc_1700-2005.dat\n");  
			exit(1); 
		}else{
            Flag_FOPEN[36] ++;
        }
		if( (fp_s[37]=fopen("./data/luc_eos2_gflsp_1700-2005.dat","rt"))==NULL ){  
			printf("No luc_eos2_gflsp_1700-2005.dat\n");  
			exit(1); 
		}else{
            Flag_FOPEN[37] ++;
        }
		if( (fp_s[38]=fopen("./data/luc_eos2_gfsh1_1700-2005.dat","rt"))==NULL ){  
			printf("No luc_eos2_gfsh1_1700-2005.dat\n");  
			exit(1); 
		}else{
            Flag_FOPEN[38] ++;
        }
		if( (fp_s[39]=fopen("./data/luc_eos2_gfsh2_1700-2005.dat","rt"))==NULL ){  
			printf("No luc_eos2_gfsh2_1700-2005.dat\n");  
			exit(1); 
		}else{
            Flag_FOPEN[39] ++;
        }
		if( (fp_s[40]=fopen("./data/luc_eos2_gfsh3_1700-2005.dat","rt"))==NULL ){  
			printf("No luc_eos2_gfsh3_1700-2005.dat\n");  
			exit(1); 
		}else{
            Flag_FOPEN[40] ++;
        }
		if( (fp_s[41]=fopen("./data/luc_eos2_gflvc_1700-2005.dat","rt"))==NULL ){  
			printf("No luc_eos2_gflvc_1700-2005.dat\n");  
			exit(1); 
		}else{
            Flag_FOPEN[41] ++;
        }
		if( (fp_s[42]=fopen("./data/luc_eos2_gflvp_1700-2005.dat","rt"))==NULL ){  
			printf("No luc_eos2_gflvp_1700-2005.dat\n");  
			exit(1); 
		}else{
            Flag_FOPEN[42] ++;
        }
		if( (fp_s[43]=fopen("./data/luc_eos2_gfvh1_1700-2005.dat","rt"))==NULL ){  
			printf("No luc_eos2_gfvh1_1700-2005.dat\n");  
			exit(1); 
		}else{
            Flag_FOPEN[43] ++;
        }
		if( (fp_s[44]=fopen("./data/luc_eos2_gfvh2_1700-2005.dat","rt"))==NULL ){  
			printf("No luc_eos2_gfvh2_1700-2005.dat\n");  
			exit(1); 
		}else{
            Flag_FOPEN[44] ++;
        }
	}
    
    /* U.NH harmonized historical land-use data (2013/12/19: by A.Ito) */
	/* URL  http://luh.unh.edu/ */
	if(LANDUSE==9 || LANDUSE==10 || LANDUSE==11 || LANDUSE==12 || LANDUSE==13
        || LANDUSE==14 || LANDUSE==15 || LANDUSE==16 || LANDUSE==17
        || LANDUSE==18|| LANDUSE == 19 || LANDUSE == 20 || LANDUSE == 21
        || LANDUSE == 22 || LANDUSE == 23 || LANDUSE == 29){
        if( (fp_s[26]=fopen("./data/luh_gcrop_1500_2005.txt","rt"))==NULL ){  
			printf("No luh_gcrop_1500_2005.txt\n");  
			exit(1); 
		}else{
            Flag_FOPEN[26] ++;
        }
		if( (fp_s[27]=fopen("./data/luh_gpast_1500_2005.txt","rt"))==NULL ){  
			printf("No luh_gpast_1500_2005.txt\n");  
			exit(1); 
		}else{
            Flag_FOPEN[27] ++;
        }
		if( (fp_s[28]=fopen("./data/luh_gothr_1500_2005.txt","rt"))==NULL ){  
			printf("No luh_gothr_1500_2005.txt\n");  
			exit(1); 
		}else{
            Flag_FOPEN[28] ++;
        }
		if( (fp_s[29]=fopen("./data/luh_gsecd_1500_2005.txt","rt"))==NULL ){  
			printf("No luh_gsecd_1500_2005.txt\n");  
			exit(1); 
		}else{
            Flag_FOPEN[29] ++;
        }
		if( (fp_s[30]=fopen("./data/luh_gssma_1500_2005.txt","rt"))==NULL ){  
			printf("No luh_gssma_1500_2005.txt\n");  
			exit(1); 
		}else{
            Flag_FOPEN[30] ++;
        }
		if( (fp_s[31]=fopen("./data/luh_gssmb_1500_2005.txt","rt"))==NULL ){  
			printf("No luh_gssmb_1500_2005.txt\n");  
			exit(1); 
		}else{
            Flag_FOPEN[31] ++;
        }
		if( (fp_s[32]=fopen("./data/luh_gflcp_1500_2005.txt","rt"))==NULL ){  
			printf("No luh_gflcp_1500_2005.txt\n");  
			exit(1); 
		}else{
            Flag_FOPEN[32] ++;
        }
		if( (fp_s[33]=fopen("./data/luh_gflcs_1500_2005.txt","rt"))==NULL ){  
			printf("No luh_gflcs_1500_2005.txt\n");  
			exit(1); 
		}else{
            Flag_FOPEN[33] ++;
        }
		if( (fp_s[34]=fopen("./data/luh_gflpc_1500_2005.txt","rt"))==NULL ){  
			printf("No luh_gflpc_1500_2005.txt\n");  
			exit(1); 
		}else{
            Flag_FOPEN[34] ++;
        }
		if( (fp_s[35]=fopen("./data/luh_gflps_1500_2005.txt","rt"))==NULL ){  
			printf("No luh_gflps_1500_2005.txt\n");  
			exit(1); 
		}else{
            Flag_FOPEN[35] ++;
        }
		if( (fp_s[36]=fopen("./data/luh_gflsc_1500_2005.txt","rt"))==NULL ){  
			printf("No luh_gflsc_1500_2005.txt\n");  
			exit(1); 
		}else{
            Flag_FOPEN[36] ++;
        }
		if( (fp_s[37]=fopen("./data/luh_gflsp_1500_2005.txt","rt"))==NULL ){  
			printf("No luh_gflsp_1500_2005.txt\n");  
			exit(1); 
		}else{
            Flag_FOPEN[37] ++;
        }
		if( (fp_s[38]=fopen("./data/luh_gfsh1_1500_2005.txt","rt"))==NULL ){  
			printf("No luh_gfsh1_1500_2005.txt\n");  
			exit(1); 
		}else{
            Flag_FOPEN[38] ++;
        }
		if( (fp_s[39]=fopen("./data/luh_gfsh2_1500_2005.txt","rt"))==NULL ){  
			printf("No luh_gfsh2_1500_2005.txt\n");  
			exit(1); 
		}else{
            Flag_FOPEN[39] ++;
        }
		if( (fp_s[40]=fopen("./data/luh_gfsh3_1500_2005.txt","rt"))==NULL ){  
			printf("No luh_gfsh3_1500_2005.txt\n");  
			exit(1); 
		}else{
            Flag_FOPEN[40] ++;
        }
		if( (fp_s[41]=fopen("./data/luh_gflvc_1500_2005.txt","rt"))==NULL ){  
			printf("No luh_gflvc_1500_2005.txt\n");  
			exit(1); 
		}else{
            Flag_FOPEN[41] ++;
        }
		if( (fp_s[42]=fopen("./data/luh_gflvp_1500_2005.txt","rt"))==NULL ){  
			printf("No luh_gflvp_1500_2005.txt\n");  
			exit(1); 
		}else{
            Flag_FOPEN[42] ++;
        }
		if( (fp_s[43]=fopen("./data/luh_gfvh1_1500_2005.txt","rt"))==NULL ){  
			printf("No luh_gfvh1_1500_2005.txt\n");  
			exit(1); 
		}else{
            Flag_FOPEN[43] ++;
        }
		if( (fp_s[44]=fopen("./data/luh_gfvh2_1500_2005.txt","rt"))==NULL ){  
			printf("No luh_gfvh2_1500_2005.txt\n");  
			exit(1); 
		}else{
            Flag_FOPEN[44] ++;
        }
    }
    
    /* land-use: projection *************************************/
    if(LANDUSE == 9 || LANDUSE == 10 || LANDUSE == 18 || LANDUSE == 19 ||
                LANDUSE == 20 || LANDUSE == 21 || LANDUSE == 22 || LANDUSE == 23){
        if( (fp_s[59]=fopen("./data/luh_gcrop_rcp45_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gcrop_rcp45_2005_2100.txt\n");  
            exit(1); 
        }else{
            Flag_FOPEN[59] ++;
        }
        if( (fp_s[60]=fopen("./data/luh_gpast_rcp45_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gpast_rcp45_2005_2100.txt\n");  
            exit(1); 
        }else{
            Flag_FOPEN[60] ++;
        }
        if( (fp_s[61]=fopen("./data/luh_gothr_rcp45_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gothr_rcp45_2005_2100.txt\n");  
            exit(1); 
        }
        if( (fp_s[62]=fopen("./data/luh_gsecd_rcp45_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gsecd_rcp45_2005_2100.txt\n");  
            exit(1); 
        }
        if( (fp_s[63]=fopen("./data/luh_gssma_rcp45_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gssma_rcp45_2005_2100.txt\n");  
            exit(1); 
        }
        if( (fp_s[64]=fopen("./data/luh_gssmb_rcp45_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gssmb_rcp45_2005_2100.txt\n");  
            exit(1); 
        }
        if( (fp_s[65]=fopen("./data/luh_gflcp_rcp45_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gflcp_rcp45_2005_2100.txt\n");  
            exit(1); 
        }
        if( (fp_s[66]=fopen("./data/luh_gflcs_rcp45_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gflcs_rcp45_2005_2100.txt\n");  
            exit(1); 
        }
        if( (fp_s[67]=fopen("./data/luh_gflpc_rcp45_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gflpc_rcp45_2005_2100.txt\n");  
            exit(1); 
        }
        if( (fp_s[68]=fopen("./data/luh_gflps_rcp45_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gflps_rcp45_2005_2100.txt\n");  
            exit(1); 
        }
        if( (fp_s[69]=fopen("./data/luh_gflsc_rcp45_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gflsc_rcp45_2005_2100.txt\n");  
            exit(1); 
        }
        if( (fp_s[70]=fopen("./data/luh_gflsp_rcp45_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gflsp_rcp45_2005_2100.txt\n");  
            exit(1); 
        }
        if( (fp_s[71]=fopen("./data/luh_gfsh1_rcp45_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gfsh1_rcp45_2005_2100.txt\n");  
            exit(1); 
        }
        if( (fp_s[72]=fopen("./data/luh_gfsh2_rcp45_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gfsh2_rcp45_2005_2100.txt\n");  
            exit(1); 
        }
        if( (fp_s[73]=fopen("./data/luh_gfsh3_rcp45_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gfsh3_rcp45_2005_2100.txt\n");  
            exit(1); 
        }
        if( (fp_s[74]=fopen("./data/luh_gflvc_rcp45_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gflvc_rcp45_2005_2100.txt\n");  
            exit(1); 
        }
        if( (fp_s[75]=fopen("./data/luh_gflvp_rcp45_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gflvp_rcp45_2005_2100.txt\n");  
            exit(1); 
        }
        if( (fp_s[76]=fopen("./data/luh_gfvh1_rcp45_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gfvh1_rcp45_2005_2100.txt\n");  
            exit(1); 
        }
        if( (fp_s[77]=fopen("./data/luh_gfvh2_rcp45_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gfvh2_rcp45_2005_2100.txt\n");  
            exit(1); 
        }
    }else if(LANDUSE == 11 || LANDUSE == 17){
        if( (fp_s[59]=fopen("./data/luh_gcrop_rcp26_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gcrop_rcp26_2005_2100.txt\n");  
            exit(1); 
        }
        if( (fp_s[60]=fopen("./data/luh_gpast_rcp26_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gpast_rcp26_2005_2100.txt\n");  
            exit(1); 
        }
        if( (fp_s[61]=fopen("./data/luh_gothr_rcp26_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gothr_rcp26_2005_2100.txt\n");  
            exit(1); 
        }
        if( (fp_s[62]=fopen("./data/luh_gsecd_rcp26_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gsecd_rcp26_2005_2100.txt\n");  
            exit(1); 
        }
        if( (fp_s[63]=fopen("./data/luh_gssma_rcp26_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gssma_rcp26_2005_2100.txt\n");  
            exit(1); 
        }
        if( (fp_s[64]=fopen("./data/luh_gssmb_rcp26_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gssmb_rcp26_2005_2100.txt\n");  
            exit(1); 
        }
        if( (fp_s[65]=fopen("./data/luh_gflcp_rcp26_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gflcp_rcp26_2005_2100.txt\n");  
            exit(1); 
        }
        if( (fp_s[66]=fopen("./data/luh_gflcs_rcp26_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gflcs_rcp26_2005_2100.txt\n");  
            exit(1); 
        }
        if( (fp_s[67]=fopen("./data/luh_gflpc_rcp26_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gflpc_rcp26_2005_2100.txt\n");  
            exit(1); 
        }
        if( (fp_s[68]=fopen("./data/luh_gflps_rcp26_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gflps_rcp26_2005_2100.txt\n");  
            exit(1); 
        }
        if( (fp_s[69]=fopen("./data/luh_gflsc_rcp26_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gflsc_rcp26_2005_2100.txt\n");  
            exit(1); 
        }
        if( (fp_s[70]=fopen("./data/luh_gflsp_rcp26_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gflsp_rcp26_2005_2100.txt\n");  
            exit(1); 
        }
        if( (fp_s[71]=fopen("./data/luh_gfsh1_rcp26_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gfsh1_rcp26_2005_2100.txt\n");  
            exit(1); 
        }
        if( (fp_s[72]=fopen("./data/luh_gfsh2_rcp26_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gfsh2_rcp26_2005_2100.txt\n");  
            exit(1); 
        }
        if( (fp_s[73]=fopen("./data/luh_gfsh3_rcp26_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gfsh3_rcp26_2005_2100.txt\n");  
            exit(1); 
        }
        if( (fp_s[74]=fopen("./data/luh_gflvc_rcp26_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gflvc_rcp26_2005_2100.txt\n");  
            exit(1); 
        }
        if( (fp_s[75]=fopen("./data/luh_gflvp_rcp26_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gflvp_rcp26_2005_2100.txt\n");  
            exit(1); 
        }
        if( (fp_s[76]=fopen("./data/luh_gfvh1_rcp26_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gfvh1_rcp26_2005_2100.txt\n");  
            exit(1); 
        }
        if( (fp_s[77]=fopen("./data/luh_gfvh2_rcp26_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gfvh2_rcp26_2005_2100.txt\n");  
            exit(1); 
        }
    }else if(LANDUSE == 12){
        if( (fp_s[59]=fopen("./data/luh_gcrop_rcp60_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gcrop_rcp60_2005_2100.txt\n");  
            exit(1); 
        }
        if( (fp_s[60]=fopen("./data/luh_gpast_rcp60_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gpast_rcp60_2005_2100.txt\n");  
            exit(1); 
        }
        if( (fp_s[61]=fopen("./data/luh_gothr_rcp60_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gothr_rcp60_2005_2100.txt\n");  
            exit(1); 
        }
        if( (fp_s[62]=fopen("./data/luh_gsecd_rcp60_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gsecd_rcp60_2005_2100.txt\n");  
            exit(1); 
        }
        if( (fp_s[63]=fopen("./data/luh_gssma_rcp60_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gssma_rcp60_2005_2100.txt\n");  
            exit(1); 
        }
        if( (fp_s[64]=fopen("./data/luh_gssmb_rcp60_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gssmb_rcp60_2005_2100.txt\n");  
            exit(1); 
        }
        if( (fp_s[65]=fopen("./data/luh_gflcp_rcp60_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gflcp_rcp60_2005_2100.txt\n");  
            exit(1); 
        }
        if( (fp_s[66]=fopen("./data/luh_gflcs_rcp60_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gflcs_rcp60_2005_2100.txt\n");  
            exit(1); 
        }
        if( (fp_s[67]=fopen("./data/luh_gflpc_rcp60_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gflpc_rcp60_2005_2100.txt\n");  
            exit(1); 
        }
        if( (fp_s[68]=fopen("./data/luh_gflps_rcp60_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gflps_rcp60_2005_2100.txt\n");  
            exit(1); 
        }
        if( (fp_s[69]=fopen("./data/luh_gflsc_rcp60_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gflsc_rcp60_2005_2100.txt\n");  
            exit(1); 
        }
        if( (fp_s[70]=fopen("./data/luh_gflsp_rcp60_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gflsp_rcp60_2005_2100.txt\n");  
            exit(1); 
        }
        if( (fp_s[71]=fopen("./data/luh_gfsh1_rcp60_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gfsh1_rcp60_2005_2100.txt\n");  
            exit(1); 
        }
        if( (fp_s[72]=fopen("./data/luh_gfsh2_rcp60_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gfsh2_rcp60_2005_2100.txt\n");  
            exit(1); 
        }
        if( (fp_s[73]=fopen("./data/luh_gfsh3_rcp60_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gfsh3_rcp60_2005_2100.txt\n");  
            exit(1); 
        }
        if( (fp_s[74]=fopen("./data/luh_gflvc_rcp60_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gflvc_rcp60_2005_2100.txt\n");  
            exit(1); 
        }
        if( (fp_s[75]=fopen("./data/luh_gflvp_rcp60_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gflvp_rcp60_2005_2100.txt\n");  
            exit(1); 
        }
        if( (fp_s[76]=fopen("./data/luh_gfvh1_rcp60_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gfvh1_rcp60_2005_2100.txt\n");  
            exit(1); 
        }
        if( (fp_s[77]=fopen("./data/luh_gfvh2_rcp60_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gfvh2_rcp60_2005_2100.txt\n");  
            exit(1); 
        }
    }else if(LANDUSE == 13){
        if( (fp_s[59]=fopen("./data/luh_gcrop_rcp85_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gcrop_rcp85_2005_2100.txt\n");  
            exit(1); 
        }
        if( (fp_s[60]=fopen("./data/luh_gpast_rcp85_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gpast_rcp85_2005_2100.txt\n");  
            exit(1); 
        }
        if( (fp_s[61]=fopen("./data/luh_gothr_rcp85_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gothr_rcp85_2005_2100.txt\n");  
            exit(1); 
        }
        if( (fp_s[62]=fopen("./data/luh_gsecd_rcp85_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gsecd_rcp85_2005_2100.txt\n");  
            exit(1); 
        }
        if( (fp_s[63]=fopen("./data/luh_gssma_rcp85_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gssma_rcp85_2005_2100.txt\n");  
            exit(1); 
        }
        if( (fp_s[64]=fopen("./data/luh_gssmb_rcp85_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gssmb_rcp85_2005_2100.txt\n");  
            exit(1); 
        }
        if( (fp_s[65]=fopen("./data/luh_gflcp_rcp85_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gflcp_rcp85_2005_2100.txt\n");  
            exit(1); 
        }
        if( (fp_s[66]=fopen("./data/luh_gflcs_rcp85_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gflcs_rcp85_2005_2100.txt\n");  
            exit(1); 
        }
        if( (fp_s[67]=fopen("./data/luh_gflpc_rcp85_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gflpc_rcp85_2005_2100.txt\n");  
            exit(1); 
        }
        if( (fp_s[68]=fopen("./data/luh_gflps_rcp85_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gflps_rcp85_2005_2100.txt\n");  
            exit(1); 
        }
        if( (fp_s[69]=fopen("./data/luh_gflsc_rcp85_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gflsc_rcp85_2005_2100.txt\n");  
            exit(1); 
        }
        if( (fp_s[70]=fopen("./data/luh_gflsp_rcp85_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gflsp_rcp85_2005_2100.txt\n");  
            exit(1); 
        }
        if( (fp_s[71]=fopen("./data/luh_gfsh1_rcp85_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gfsh1_rcp85_2005_2100.txt\n");  
            exit(1); 
        }
        if( (fp_s[72]=fopen("./data/luh_gfsh2_rcp85_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gfsh2_rcp85_2005_2100.txt\n");  
            exit(1); 
        }
        if( (fp_s[73]=fopen("./data/luh_gfsh3_rcp85_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gfsh3_rcp85_2005_2100.txt\n");  
            exit(1); 
        }
        if( (fp_s[74]=fopen("./data/luh_gflvc_rcp85_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gflvc_rcp85_2005_2100.txt\n");  
            exit(1); 
        }
        if( (fp_s[75]=fopen("./data/luh_gflvp_rcp85_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gflvp_rcp85_2005_2100.txt\n");  
            exit(1); 
        }
        if( (fp_s[76]=fopen("./data/luh_gfvh1_rcp85_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gfvh1_rcp85_2005_2100.txt\n");  
            exit(1); 
        }
        if( (fp_s[77]=fopen("./data/luh_gfvh2_rcp85_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gfvh2_rcp85_2005_2100.txt\n");  
            exit(1); 
        }
    }
    
    /* ICARUS SSP by T.Kinoshita: 2014/09/04 ******************/
    /* SSP1 */
    if(LANDUSE == 14){
        if( (fp_s[59]=fopen("./data/luc_crop_ssp1_2001-2100.txt","rt"))==NULL ){
            printf("No luc_crop_ssp1_2001-2100.txt\n");
            exit(1); 
        }
        if( (fp_s[60]=fopen("./data/luc_pasture_ssp1_2001-2100.txt","rt"))==NULL ){
            printf("No luc_pasture_ssp1_2001-2100.txt\n");
            exit(1); 
        }
        if( (fp_s[61]=fopen("./data/luc_other_ssp1_2001-2100.txt","rt"))==NULL ){
            printf("No luc_other_ssp1_2001-2100.txt\n");
            exit(1); 
        }
        if( (fp_s[62]=fopen("./data/luc_mf_ssp1_2001-2100.txt","rt"))==NULL ){
            printf("No luc_mf_ssp1_2001-2100.txt\n");
            exit(1); 
        }
        if( (fp_s[63]=fopen("./data/luc_grass_ssp1_2001-2100.txt","rt"))==NULL ){
            printf("No luc_grass_ssp1_2001-2100.txt\n");
            exit(1); 
        }
        if( (fp_s[64]=fopen("./data/luc_umf_ssp1_2001-2100.txt","rt"))==NULL ){
            printf("No luc_umf_ssp1_2001-2100.txt\n");
            exit(1); 
        }
        if( (fp_s[65]=fopen("./data/luc_c2p_ssp1_2001-2100.txt","rt"))==NULL ){
            printf("No luc_c2p_ssp1_2001-2100.txt\n");
            exit(1); 
        }
        if( (fp_s[66]=fopen("./data/luc_c2s_ssp1_2001-2100.txt","rt"))==NULL ){
            printf("No luc_c2s_ssp1_2001-2100.txt\n");
            exit(1); 
        }
        if( (fp_s[67]=fopen("./data/luc_c2u_ssp1_2001-2100.txt","rt"))==NULL ){
            printf("No luc_c2u_ssp1_2001-2100.txt\n");
            exit(1); 
        }
        if( (fp_s[68]=fopen("./data/luc_p2c_ssp1_2001-2100.txt","rt"))==NULL ){
            printf("No luc_p2c_ssp1_2001-2100.txt\n");
            exit(1); 
        }
        if( (fp_s[69]=fopen("./data/luc_p2s_ssp1_2001-2100.txt","rt"))==NULL ){
            printf("No luc_p2s_ssp1_2001-2100.txt\n");
            exit(1); 
        }
        if( (fp_s[70]=fopen("./data/luc_p2u_ssp1_2001-2100.txt","rt"))==NULL ){
            printf("No luc_p2u_ssp1_2001-2100.txt\n");
            exit(1); 
        }
        if( (fp_s[71]=fopen("./data/luc_s2c_ssp1_2001-2100.txt","rt"))==NULL ){
            printf("No luc_s2c_ssp1_2001-2100.txt\n");
            exit(1); 
        }
        if( (fp_s[72]=fopen("./data/luc_s2p_ssp1_2001-2100.txt","rt"))==NULL ){
            printf("No luc_s2p_ssp1_2001-2100.txt\n");
            exit(1); 
        }
        if( (fp_s[73]=fopen("./data/luc_s2u_ssp1_2001-2100.txt","rt"))==NULL ){
            printf("No luc_s2u_ssp1_2001-2100.txt\n");
            exit(1); 
        }
        if( (fp_s[74]=fopen("./data/luc_v2c_ssp1_2001-2100.txt","rt"))==NULL ){
            printf("No luc_v2c_ssp1_2001-2100.txt\n");
            exit(1); 
        }
        if( (fp_s[75]=fopen("./data/luc_v2p_ssp1_2001-2100.txt","rt"))==NULL ){
            printf("No luc_v2p_ssp1_2001-2100.txt\n");
            exit(1); 
        }
        if( (fp_s[76]=fopen("./data/luc_v2s_ssp1_2001-2100.txt","rt"))==NULL ){
            printf("No luc_v2s_ssp1_2001-2100.txt\n");
            exit(1); 
        }
        if( (fp_s[77]=fopen("./data/luc_v2u_ssp1_2001-2100.txt","rt"))==NULL ){
            printf("No luc_v2u_ssp1_2001-2100.txt\n");
            exit(1); 
        }
    }
    
    /* SSP2 */
    if(LANDUSE == 15){
        if( (fp_s[59]=fopen("./data/luc_crop_ssp2_2001-2100.txt","rt"))==NULL ){
            printf("No luc_crop_ssp2_2001-2100.txt\n");
            exit(1); 
        }
        if( (fp_s[60]=fopen("./data/luc_pasture_ssp2_2001-2100.txt","rt"))==NULL ){
            printf("No luc_pasture_ssp2_2001-2100.txt\n");
            exit(1); 
        }
        if( (fp_s[61]=fopen("./data/luc_other_ssp2_2001-2100.txt","rt"))==NULL ){
            printf("No luc_other_ssp2_2001-2100.txt\n");
            exit(1); 
        }
        if( (fp_s[62]=fopen("./data/luc_mf_ssp2_2001-2100.txt","rt"))==NULL ){
            printf("No luc_mf_ssp2_2001-2100.txt\n");
            exit(1); 
        }
        if( (fp_s[63]=fopen("./data/luc_grass_ssp2_2001-2100.txt","rt"))==NULL ){
            printf("No luc_grass_ssp2_2001-2100.txt\n");
            exit(1); 
        }
        if( (fp_s[64]=fopen("./data/luc_umf_ssp2_2001-2100.txt","rt"))==NULL ){
            printf("No luc_umf_ssp2_2001-2100.txt\n");
            exit(1); 
        }
        if( (fp_s[65]=fopen("./data/luc_c2p_ssp2_2001-2100.txt","rt"))==NULL ){
            printf("No luc_c2p_ssp2_2001-2100.txt\n");
            exit(1); 
        }
        if( (fp_s[66]=fopen("./data/luc_c2s_ssp2_2001-2100.txt","rt"))==NULL ){
            printf("No luc_c2s_ssp2_2001-2100.txt\n");
            exit(1); 
        }
        if( (fp_s[67]=fopen("./data/luc_c2u_ssp2_2001-2100.txt","rt"))==NULL ){
            printf("No luc_c2u_ssp2_2001-2100.txt\n");
            exit(1); 
        }
        if( (fp_s[68]=fopen("./data/luc_p2c_ssp2_2001-2100.txt","rt"))==NULL ){
            printf("No luc_p2c_ssp2_2001-2100.txt\n");
            exit(1); 
        }
        if( (fp_s[69]=fopen("./data/luc_p2s_ssp2_2001-2100.txt","rt"))==NULL ){
            printf("No luc_p2s_ssp2_2001-2100.txt\n");
            exit(1); 
        }
        if( (fp_s[70]=fopen("./data/luc_p2u_ssp2_2001-2100.txt","rt"))==NULL ){
            printf("No luc_p2u_ssp2_2001-2100.txt\n");
            exit(1); 
        }
        if( (fp_s[71]=fopen("./data/luc_s2c_ssp2_2001-2100.txt","rt"))==NULL ){
            printf("No luc_s2c_ssp2_2001-2100.txt\n");
            exit(1); 
        }
        if( (fp_s[72]=fopen("./data/luc_s2p_ssp2_2001-2100.txt","rt"))==NULL ){
            printf("No luc_s2p_ssp2_2001-2100.txt\n");
            exit(1); 
        }
        if( (fp_s[73]=fopen("./data/luc_s2u_ssp2_2001-2100.txt","rt"))==NULL ){
            printf("No luc_s2u_ssp2_2001-2100.txt\n");
            exit(1); 
        }
        if( (fp_s[74]=fopen("./data/luc_v2c_ssp2_2001-2100.txt","rt"))==NULL ){
            printf("No luc_v2c_ssp2_2001-2100.txt\n");
            exit(1); 
        }
        if( (fp_s[75]=fopen("./data/luc_v2p_ssp2_2001-2100.txt","rt"))==NULL ){
            printf("No luc_v2p_ssp2_2001-2100.txt\n");
            exit(1); 
        }
        if( (fp_s[76]=fopen("./data/luc_v2s_ssp2_2001-2100.txt","rt"))==NULL ){
            printf("No luc_v2s_ssp2_2001-2100.txt\n");
            exit(1); 
        }
        if( (fp_s[77]=fopen("./data/luc_v2u_ssp2_2001-2100.txt","rt"))==NULL ){
            printf("No luc_v2u_ssp2_2001-2100.txt\n");
            exit(1); 
        }
    }
    
    /* SSP3 */
    if(LANDUSE == 16){
        if( (fp_s[59]=fopen("./data/luc_crop_ssp3_2001-2100.txt","rt"))==NULL ){
            printf("No luc_crop_ssp3_2001-2100.txt\n");
            exit(1); 
        }
        if( (fp_s[60]=fopen("./data/luc_pasture_ssp3_2001-2100.txt","rt"))==NULL ){
            printf("No luc_pasture_ssp3_2001-2100.txt\n");
            exit(1); 
        }
        if( (fp_s[61]=fopen("./data/luc_other_ssp3_2001-2100.txt","rt"))==NULL ){
            printf("No luc_other_ssp3_2001-2100.txt\n");
            exit(1); 
        }
        if( (fp_s[62]=fopen("./data/luc_mf_ssp3_2001-2100.txt","rt"))==NULL ){
            printf("No luc_mf_ssp3_2001-2100.txt\n");
            exit(1); 
        }
        if( (fp_s[63]=fopen("./data/luc_grass_ssp3_2001-2100.txt","rt"))==NULL ){
            printf("No luc_grass_ssp3_2001-2100.txt\n");
            exit(1); 
        }
        if( (fp_s[64]=fopen("./data/luc_umf_ssp3_2001-2100.txt","rt"))==NULL ){
            printf("No luc_umf_ssp3_2001-2100.txt\n");
            exit(1); 
        }
        if( (fp_s[65]=fopen("./data/luc_c2p_ssp3_2001-2100.txt","rt"))==NULL ){
            printf("No luc_c2p_ssp3_2001-2100.txt\n");
            exit(1); 
        }
        if( (fp_s[66]=fopen("./data/luc_c2s_ssp3_2001-2100.txt","rt"))==NULL ){
            printf("No luc_c2s_ssp3_2001-2100.txt\n");
            exit(1); 
        }
        if( (fp_s[67]=fopen("./data/luc_c2u_ssp3_2001-2100.txt","rt"))==NULL ){
            printf("No luc_c2u_ssp3_2001-2100.txt\n");
            exit(1); 
        }
        if( (fp_s[68]=fopen("./data/luc_p2c_ssp3_2001-2100.txt","rt"))==NULL ){
            printf("No luc_p2c_ssp3_2001-2100.txt\n");
            exit(1); 
        }
        if( (fp_s[69]=fopen("./data/luc_p2s_ssp3_2001-2100.txt","rt"))==NULL ){
            printf("No luc_p2s_ssp3_2001-2100.txt\n");
            exit(1); 
        }
        if( (fp_s[70]=fopen("./data/luc_p2u_ssp3_2001-2100.txt","rt"))==NULL ){
            printf("No luc_p2u_ssp3_2001-2100.txt\n");
            exit(1); 
        }
        if( (fp_s[71]=fopen("./data/luc_s2c_ssp3_2001-2100.txt","rt"))==NULL ){
            printf("No luc_s2c_ssp3_2001-2100.txt\n");
            exit(1); 
        }
        if( (fp_s[72]=fopen("./data/luc_s2p_ssp3_2001-2100.txt","rt"))==NULL ){
            printf("No luc_s2p_ssp3_2001-2100.txt\n");
            exit(1); 
        }
        if( (fp_s[73]=fopen("./data/luc_s2u_ssp3_2001-2100.txt","rt"))==NULL ){
            printf("No luc_s2u_ssp3_2001-2100.txt\n");
            exit(1); 
        }
        if( (fp_s[74] = fopen("./data/luc_v2c_ssp3_2001-2100.txt","rt"))==NULL ){
            printf("No luc_v2c_ssp3_2001-2100.txt\n");
            exit(1); 
        }
        if( (fp_s[75] = fopen("./data/luc_v2p_ssp3_2001-2100.txt","rt"))==NULL ){
            printf("No luc_v2p_ssp3_2001-2100.txt\n");
            exit(1); 
        }
        if( (fp_s[76] = fopen("./data/luc_v2s_ssp3_2001-2100.txt","rt"))==NULL ){
            printf("No luc_v2s_ssp3_2001-2100.txt\n");
            exit(1); 
        }
        if( (fp_s[77] = fopen("./data/luc_v2u_ssp3_2001-2100.txt","rt"))==NULL ){
            printf("No luc_v2u_ssp3_2001-2100.txt\n");
            exit(1); 
        }
    }
    
    /* TELUMO LUC (2019/02/21 by A.Ito) */
    if(LANDUSE == 29){
        switch(SCENARIO_ID){
            case 2001: strcpy(telumoid,"./h26/"); break;
            case 2011: strcpy(telumoid,"./i26/"); break;
            case 2021: strcpy(telumoid,"./g26/"); break;
            case 2031: strcpy(telumoid,"./m26/"); break;
            case 2041: strcpy(telumoid,"./n26/"); break;
            case 2002: strcpy(telumoid,"./h85/"); break;
            case 2012: strcpy(telumoid,"./i85/"); break;
            case 2022: strcpy(telumoid,"./g85/"); break;
            case 2032: strcpy(telumoid,"./m85/"); break;
            case 2042: strcpy(telumoid,"./n85/"); break;
            case 2003: strcpy(telumoid,"./h45/"); break;
            case 2013: strcpy(telumoid,"./i45/"); break;
            case 2023: strcpy(telumoid,"./g45/"); break;
            case 2033: strcpy(telumoid,"./m45/"); break;
            case 2043: strcpy(telumoid,"./n45/"); break;
            case 2004: strcpy(telumoid,"./h85/"); break;
            case 2014: strcpy(telumoid,"./i85/"); break;
            case 2024: strcpy(telumoid,"./g85/"); break;
            case 2034: strcpy(telumoid,"./m85/"); break;
            case 2044: strcpy(telumoid,"./n85/"); break;
            default:
                printf("!!!! BAD TELUMO LUC ID\n"); exit(1); break;
        }
        
        strcpy(filename,telumoid); strcat(filename,"luc_crop_telmo_2005-2100.txt");
        if( (fp_s[59]=fopen(filename,"rt"))==NULL ){
            printf("No luc_crop_telmo_2005-2100.txt\n");
            exit(1);
        }else{
            Flag_FOPEN[59] ++;
        }
        
        strcpy(filename,telumoid); strcat(filename,"luc_pasture_telmo_2005-2100.txt");
        if( (fp_s[60]=fopen(filename,"rt"))==NULL ){
            printf("No luc_pasture_telmo_2005-2100.txt\n");
            exit(1);
        }else{
            Flag_FOPEN[60] ++;
        }
        
        strcpy(filename,telumoid); strcat(filename,"luc_other_telmo_2005-2100.txt");
        if( (fp_s[61]=fopen(filename,"rt"))==NULL ){
            printf("No luc_other_telmo_2005-2100.txt\n");
            exit(1);
        }else{
            Flag_FOPEN[61] ++;
        }

        strcpy(filename,telumoid); strcat(filename,"luc_mf_telmo_2005-2100.txt");
        if( (fp_s[62]=fopen(filename,"rt"))==NULL ){
            printf("No luc_mf_telmo_2005-2100.txt\n");
            exit(1);
        }else{
            Flag_FOPEN[62] ++;
        }

        strcpy(filename,telumoid); strcat(filename,"luc_grass_telmo_2005-2100.txt");
        if( (fp_s[63]=fopen(filename,"rt"))==NULL ){
            printf("No luc_grass_telmo_2005-2100.txt\n");
            exit(1);
        }else{
            Flag_FOPEN[63] ++;
        }

        strcpy(filename,telumoid); strcat(filename,"luc_umf_telmo_2005-2100.txt");
        if( (fp_s[64]=fopen(filename,"rt"))==NULL ){
            printf("No luc_umf_telmo_2005-2100.txt\n");
            exit(1);
        }else{
            Flag_FOPEN[64] ++;
        }

        strcpy(filename,telumoid); strcat(filename,"luc_c2p_telmo_2005-2100.txt");
        if( (fp_s[65]=fopen(filename,"rt"))==NULL ){
            printf("No luc_c2p_telmo_2005-2100.txt\n");
            exit(1);
        }else{
            Flag_FOPEN[65] ++;
        }

        strcpy(filename,telumoid); strcat(filename,"luc_c2s_telmo_2005-2100.txt");
        if( (fp_s[66]=fopen(filename,"rt"))==NULL ){
            printf("No luc_c2s_telmo_2005-2100.txt\n");
            exit(1);
        }else{
            Flag_FOPEN[66] ++;
        }

        strcpy(filename,telumoid); strcat(filename,"luc_c2u_telmo_2005-2100.txt");
        if( (fp_s[67]=fopen(filename,"rt"))==NULL ){
            printf("No luc_c2u_telmo_2005-2100.txt\n");
            exit(1);
        }else{
            Flag_FOPEN[67] ++;
        }

        strcpy(filename,telumoid); strcat(filename,"luc_p2c_telmo_2005-2100.txt");
        if( (fp_s[68]=fopen(filename,"rt"))==NULL ){
            printf("No luc_p2c_telmo_2005-2100.txt\n");
            exit(1);
        }else{
            Flag_FOPEN[68] ++;
        }

        strcpy(filename,telumoid); strcat(filename,"luc_p2s_telmo_2005-2100.txt");
        if( (fp_s[69]=fopen(filename,"rt"))==NULL ){
            printf("No luc_p2s_telmo_2005-2100.txt\n");
            exit(1);
        }else{
            Flag_FOPEN[69] ++;
        }

        strcpy(filename,telumoid); strcat(filename,"luc_p2u_telmo_2005-2100.txt");
        if( (fp_s[70]=fopen(filename,"rt"))==NULL ){
            printf("No luc_p2u_telmo_2005-2100.txt\n");
            exit(1);
        }else{
            Flag_FOPEN[70] ++;
        }

        strcpy(filename,telumoid); strcat(filename,"luc_s2c_telmo_2005-2100.txt");
        if( (fp_s[71]=fopen(filename,"rt"))==NULL ){
            printf("No luc_s2c_telmo_2005-2100.txt\n");
            exit(1);
        }else{
            Flag_FOPEN[71] ++;
        }

        strcpy(filename,telumoid); strcat(filename,"luc_s2p_telmo_2005-2100.txt");
        if( (fp_s[72]=fopen(filename,"rt"))==NULL ){
            printf("No luc_s2p_telmo_2005-2100.txt\n");
            exit(1);
        }else{
            Flag_FOPEN[72] ++;
        }
        
        strcpy(filename,telumoid); strcat(filename,"luc_s2u_telmo_2005-2100.txt");
        if( (fp_s[73]=fopen(filename,"rt"))==NULL ){
            printf("No luc_s2u_telmo_2005-2100.txt\n");
            exit(1);
        }else{
            Flag_FOPEN[73] ++;
        }

        strcpy(filename,telumoid); strcat(filename,"luc_v2c_telmo_2005-2100.txt");
        if( (fp_s[74]=fopen(filename,"rt"))==NULL ){
            printf("No luc_v2c_telmo_2005-2100.txt\n");
            exit(1);
        }else{
            Flag_FOPEN[74] ++;
        }

        strcpy(filename,telumoid); strcat(filename,"luc_v2p_telmo_2005-2100.txt");
        if( (fp_s[75]=fopen(filename,"rt"))==NULL ){
            printf("No luc_v2p_telmo_2005-2100.txt\n");
            exit(1);
        }else{
            Flag_FOPEN[75] ++;
        }

        strcpy(filename,telumoid); strcat(filename,"luc_v2s_telmo_2005-2100.txt");
        if( (fp_s[76]=fopen(filename,"rt"))==NULL ){
            printf("No luc_v2s_telmo_2005-2100.txt\n");
            exit(1);
        }else{
            Flag_FOPEN[76] ++;
        }

        strcpy(filename,telumoid); strcat(filename,"luc_v2u_telmo_2005-2100.txt");
        if( (fp_s[77]=fopen(filename,"rt"))==NULL ){
            printf("No luc_v2u_telmo_2005-2100.txt\n");
            exit(1);
        }else{
            Flag_FOPEN[77] ++;
        }
    }
    
    /* CMIP6 historical data: 1886-2015: 2018/12/21 by A.Ito */
    if(LANDUSE == 26){
        if( (fp_s[26]=fopen("./data/luh2_state_5-05_1866-2015.txt","rt"))==NULL ){
            printf("No luh2_state_5-05_1866-2015.txt\n");
            exit(1);
        }else{
            Flag_FOPEN[26] ++;
        }
        if( (fp_s[27]=fopen("./data/luh2_transition_5-05_1866-2015.txt","rt"))==NULL ){
            printf("No luh2_transition_5-05_1866-2015.txt\n");
            exit(1);
        }else{
            Flag_FOPEN[27] ++;
        }
    }
    if(LANDUSE == 27){
        if( (fp_s[26]=fopen("./data/luh2_state_5-05_1866-2015_high.txt","rt"))==NULL ){
            printf("No luh2_state_5-05_1866-2015_high.txt\n");
            exit(1);
        }else{
            Flag_FOPEN[26] ++;
        }
        if( (fp_s[27]=fopen("./data/luh2_transition_5-05_1866-2015_high.txt","rt"))==NULL ){
            printf("No luh2_transition_5-05_1866-2015_high.txt\n");
            exit(1);
        }else{
            Flag_FOPEN[27] ++;
        }
    }
    if(LANDUSE == 28){
        if( (fp_s[26]=fopen("./data/luh2_state_5-05_1866-2015_low.txt","rt"))==NULL ){
            printf("No luh2_state_5-05_1866-2015_low.txt\n");
            exit(1);
        }else{
            Flag_FOPEN[26] ++;
        }
        if( (fp_s[27]=fopen("./data/luh2_transition_5-05_1866-2015_low.txt","rt"))==NULL ){
            printf("No luh2_transition_5-05_1866-2015_low.txt\n");
            exit(1);
        }else{
            Flag_FOPEN[27] ++;
        }
    }

	/***************************************************/
	/* 0: stable */
	/* 1: gradual rise (SRES A1 or A1B) */
	/* 2: gradual rise (SRES A1FI) */
	/* 3: gradual rise (SRES A1T) */
	/* 4: gradual rise (SRES A2) */
	/* 5: gradual rise (SRES B1) */
	/* 6: gradual rise (SRES B2) */
	/* 7: uniform rise */
	if(SCENARIO_ID == 0 || (SCENARIO_ID == 6001 || SCENARIO_ID == 6002)){
		CO2S = 4;
		fp_s[23] = fopen("./data/image_a2_fcrop.dat","rt");
		fp_s[45] = fopen("./data/image_a2_fgrass.dat","rt");
	}else if(SCENARIO_ID ==1 || SCENARIO_ID == 41 || SCENARIO_ID == 61){
		CO2S = 1;
		fp_s[23]=fopen("./data/image_a1b_fcrop.dat","rt");
		fp_s[45]=fopen("./data/image_a1b_fgrass.dat","rt");
	}else if(SCENARIO_ID == 2){
		CO2S = 2;
		fp_s[23]=fopen("./data/image_a1b_fcrop.dat","rt");
		fp_s[45]=fopen("./data/image_a1b_fgrass.dat","rt");
	}else if(SCENARIO_ID == 3){
		CO2S = 3;
		fp_s[23]=fopen("./data/image_a1b_fcrop.dat","rt");
		fp_s[45]=fopen("./data/image_a1b_fgrass.dat","rt");
	}else if(SCENARIO_ID==4 || SCENARIO_ID==11 || SCENARIO_ID==13 || SCENARIO_ID==14 || SCENARIO_ID==15 ||
            SCENARIO_ID==21 || SCENARIO_ID==22 || SCENARIO_ID==23 || SCENARIO_ID==31 ||
			SCENARIO_ID==42 || SCENARIO_ID==51 || SCENARIO_ID==62 || SCENARIO_ID==71 ){
		CO2S = 4;
		fp_s[23]=fopen("./data/image_a2_fcrop.dat","rt");
		fp_s[45]=fopen("./data/image_a2_fgrass.dat","rt");
	}else if(SCENARIO_ID==5 || SCENARIO_ID==43){
		CO2S = 5;
		fp_s[23]=fopen("./data/image_b1_fcrop.dat","rt");
		fp_s[45]=fopen("./data/image_b1_fgrass.dat","rt");
	}else if(SCENARIO_ID==6 || SCENARIO_ID==12 || SCENARIO_ID==16 || SCENARIO_ID==17 ||
            SCENARIO_ID==18 || SCENARIO_ID==24 || SCENARIO_ID==32 ||
			SCENARIO_ID==44 || SCENARIO_ID==52 || SCENARIO_ID==63 || SCENARIO_ID==106){
		CO2S = 6;
		fp_s[23]=fopen("./data/image_b1_fcrop.dat","rt");
		fp_s[45]=fopen("./data/image_b1_fgrass.dat","rt");
	}
	
	if(SCENARIO_ID==1000 ||SCENARIO_ID==1010 ||SCENARIO_ID==1011 ||SCENARIO_ID==1012 ||SCENARIO_ID==1050 ||
        SCENARIO_ID==1060 ||SCENARIO_ID==1070 ||
		SCENARIO_ID==1080 ||SCENARIO_ID==1090 ||SCENARIO_ID==1100 ||SCENARIO_ID==1101 ||SCENARIO_ID==1110 ||SCENARIO_ID==1120 ||
		SCENARIO_ID==1130 ||SCENARIO_ID==1131 ||SCENARIO_ID==1132 ||SCENARIO_ID==1133 ||SCENARIO_ID==1134 ||SCENARIO_ID==1150 ||
		SCENARIO_ID==1151 ||SCENARIO_ID==1152 ||SCENARIO_ID==1153 ||SCENARIO_ID==1160 ||SCENARIO_ID==1161 ||SCENARIO_ID==1162 ||
		SCENARIO_ID==1170 ||SCENARIO_ID==1180 ||SCENARIO_ID==1190 ||SCENARIO_ID==1200 ||SCENARIO_ID==1201 ||SCENARIO_ID==1202 ||
		SCENARIO_ID==1203 ||SCENARIO_ID==1204 ||SCENARIO_ID==1220 ||SCENARIO_ID==1221 ||SCENARIO_ID==1222 ||SCENARIO_ID==1223 ||
		SCENARIO_ID==1224 ||SCENARIO_ID==1225 ||SCENARIO_ID==1226 ||SCENARIO_ID==1240 ||SCENARIO_ID==1250 ||SCENARIO_ID==1251 ||
		SCENARIO_ID==1260 ||SCENARIO_ID==1261 ||SCENARIO_ID==1262 ||SCENARIO_ID==1263 ||SCENARIO_ID==1270){ 
		CO2S = 1;
		
		if((fp_s[23]=fopen("./data/image_a1b_fcrop.dat","rt"))==NULL){
			printf("NO DATA:image_a1b_fcrop.dat !!\n");
			exit(1);
		}else{
            Flag_FOPEN[23] ++;
        }
        if((fp_s[45]=fopen("./data/image_a1b_fgrass.dat","rt"))==NULL){
			printf("NO DATA:image_a1b_fgrass.dat !!\n");
			exit(1);
        }else{
            Flag_FOPEN[45] ++;
        }
	}else if(SCENARIO_ID==1013 ||SCENARIO_ID==1014 ||SCENARIO_ID==1015 || SCENARIO_ID==1051 ||
        SCENARIO_ID==1061 ||SCENARIO_ID==1071 ||SCENARIO_ID==1081 ||SCENARIO_ID==1091 ||
        SCENARIO_ID==1111 ||SCENARIO_ID==1121 ||SCENARIO_ID==1135 ||SCENARIO_ID==1136 ||
        SCENARIO_ID==1137 ||SCENARIO_ID==1138 ||SCENARIO_ID==1139 ||SCENARIO_ID==1154 ||
        SCENARIO_ID==1155 ||SCENARIO_ID==1156 ||SCENARIO_ID==1171 ||SCENARIO_ID==1181 ||
		SCENARIO_ID==1191 ||SCENARIO_ID==1205 ||SCENARIO_ID==1206 ||SCENARIO_ID==1207 ||
        SCENARIO_ID==1208 ||SCENARIO_ID==1209 ||SCENARIO_ID==1227 ||SCENARIO_ID==1228 ||
        SCENARIO_ID==1229 ||SCENARIO_ID==1230 ||SCENARIO_ID==1231 ||SCENARIO_ID==1241 ||
        SCENARIO_ID==1252 ||SCENARIO_ID==1264 ||SCENARIO_ID==1265 ||SCENARIO_ID==1266 ||
        SCENARIO_ID==1267 ||SCENARIO_ID==1271){
		CO2S = 4;

		if((fp_s[23]=fopen("./data/image_a2_fcrop.dat","rt"))==NULL){
			printf("NO DATA !!\n");
			exit(1);
		}else{
            Flag_FOPEN[23] ++;
        }
        if((fp_s[45]=fopen("./data/image_a2_fgrass.dat","rt"))==NULL){
			printf("NO DATA:image_a2_fgrass.dat !!\n");
			exit(1);
        }else{
            Flag_FOPEN[45] ++;
        }
	}else if(SCENARIO_ID==1001 ||SCENARIO_ID==1016 ||SCENARIO_ID==1017 ||SCENARIO_ID==1018 ||
        SCENARIO_ID==1052 ||SCENARIO_ID==1062 ||SCENARIO_ID==1072 ||SCENARIO_ID==1082 ||
        SCENARIO_ID==1092 ||SCENARIO_ID==1102 ||SCENARIO_ID==1103 ||SCENARIO_ID==1112 ||
        SCENARIO_ID==1122 ||SCENARIO_ID==1140 ||SCENARIO_ID==1141 ||SCENARIO_ID==1142 ||
        SCENARIO_ID==1143 ||SCENARIO_ID==1144 ||SCENARIO_ID==1157 ||SCENARIO_ID==1158 ||
        SCENARIO_ID==1159 ||SCENARIO_ID==1163 ||SCENARIO_ID==1164 ||SCENARIO_ID==1165 ||
        SCENARIO_ID==1182 ||SCENARIO_ID==1192 ||SCENARIO_ID==1210 ||SCENARIO_ID==1211 ||
		SCENARIO_ID==1212 ||SCENARIO_ID==1213 ||SCENARIO_ID==1214 ||SCENARIO_ID==1232 ||
        SCENARIO_ID==1233 ||SCENARIO_ID==1234 ||SCENARIO_ID==1235 ||SCENARIO_ID==1236 ||
        SCENARIO_ID==1237 ||SCENARIO_ID==1238 ||SCENARIO_ID==1239 ||SCENARIO_ID==1242 ||
        SCENARIO_ID==1253 ||SCENARIO_ID==1272){
		CO2S = 5;
		
		if((fp_s[23]=fopen("./data/image_b1_fcrop.dat","rt"))==NULL){
			printf("NO DATA !!\n");
			exit(1);
		}else{
            Flag_FOPEN[23] ++;
        }
        if((fp_s[45] = fopen("./data/image_b1_fgrass.dat","rt"))==NULL){
			printf("NO DATA:image_b1_fgrass.dat !!\n");
			exit(1);
        }else{
            Flag_FOPEN[45] ++;
        }
	}
    
    /* CD-LINKS: 2016/11/17 by A.Ito */
    if(LANDUSE == 19){
        if((fp_s[23] = fopen("./data/fcrop_ipsl_aim_ssp1_rcp45_v1.txt","rt"))==NULL){
            printf("NO ./data/fcrop_ipsl_aim_ssp1_rcp45_v1.txt !!\n");
            exit(1);
        }
        fp_s[45]=fopen("./data/image_b1_fgrass.dat","rt");
    }else if(LANDUSE == 20){
        if((fp_s[23]=fopen("./data/fcrop_ipsl_aim_ssp2_rcp45_v1.txt","rt"))==NULL){
            printf("NO ./data/fcrop_ipsl_aim_ssp1_rcp45_v1.txt !!\n");
            exit(1);
        }
        fp_s[45]=fopen("./data/image_b1_fgrass.dat","rt");
    }else if(LANDUSE == 21){
        if((fp_s[23]=fopen("./data/fcrop_ipsl_aim_ssp3_rcp45_v1.txt","rt"))==NULL){
            printf("NO ./data/fcrop_ipsl_aim_ssp3_rcp45_v1.txt !!\n");
            exit(1);
        }
        fp_s[45]=fopen("./data/image_b1_fgrass.dat","rt");
    }else if(LANDUSE == 22){
        if((fp_s[23]=fopen("./data/fcrop_ipsl_aim_ssp4_rcp45_v1.txt","rt"))==NULL){
            printf("NO ./data/fcrop_ipsl_aim_ssp4_rcp45_v1.txt !!\n");
            exit(1);
        }
        fp_s[45]=fopen("./data/image_b1_fgrass.dat","rt");
    }else if(LANDUSE == 23){
        if((fp_s[23]=fopen("./data/fcrop_ipsl_aim_ssp5_rcp45_v1.txt","rt"))==NULL){
            printf("NO ./data/fcrop_ipsl_aim_ssp5_rcp45_v1.txt !!\n");
            exit(1);
        }
        fp_s[45]=fopen("./data/image_b1_fgrass.dat","rt");
    }else{
        /* ISI-MIP: 2012/06/27 by A.Ito ***********************************/
        /* assumption: SRES data were used, instead of RCP data */
        if(SCENARIO_ID==2001 || SCENARIO_ID==2011 || SCENARIO_ID==2021 || SCENARIO_ID==2031 || SCENARIO_ID==2041
            || SCENARIO_ID==2005 || SCENARIO_ID==2015 || SCENARIO_ID==2025 || SCENARIO_ID==2035 || SCENARIO_ID==2045
            || SCENARIO_ID==4023){
            CO2S = 1; /* RCP2.6 */
            /**/
            if((fp_s[23]=fopen("./data/image_b1_fcrop.dat","rt"))==NULL){
                printf("NO image_b1_fcrop.dat !!\n");
                exit(1);
            }
            fp_s[45]=fopen("./data/image_b1_fgrass.dat","rt");
        }else if(SCENARIO_ID==2002 || SCENARIO_ID==2012 || SCENARIO_ID==2022 || SCENARIO_ID==2032 || SCENARIO_ID==2042
            || SCENARIO_ID==2006 || SCENARIO_ID==2016 || SCENARIO_ID==2026 || SCENARIO_ID==2036 || SCENARIO_ID==2046
            || SCENARIO_ID==4012 || SCENARIO_ID==4022){
            CO2S = 4; /* RCP8.5 */
            /**/
            if((fp_s[23] = fopen("./data/image_a2_fcrop.dat","rt"))==NULL){
                printf("NO image_a2_fcrop.dat !!\n");
                exit(1);
            }
            fp_s[45]=fopen("./data/image_a2_fgrass.dat","rt");
        }else if(SCENARIO_ID==2003 || SCENARIO_ID==2013 || SCENARIO_ID==2023 || SCENARIO_ID==2033 || SCENARIO_ID==2043
            || SCENARIO_ID==2007 || SCENARIO_ID==2017 || SCENARIO_ID==2027 || SCENARIO_ID==2037 || SCENARIO_ID==2047
            || SCENARIO_ID==4011|| SCENARIO_ID==4021){
            CO2S = 2; /* RCP4.5 */
            /**/
            if((fp_s[23]=fopen("./data/image_a1b_fcrop.dat","rt"))==NULL){
                printf("NO image_a1b_fcrop.dat !!\n");
                exit(1);
            }
            fp_s[45]=fopen("./data/image_a1b_fgrass.dat","rt");
        }else if(SCENARIO_ID==2004 || SCENARIO_ID==2014 || SCENARIO_ID==2024 || SCENARIO_ID==2034 || SCENARIO_ID==2044
            || SCENARIO_ID==2008 || SCENARIO_ID==2018 || SCENARIO_ID==2028 || SCENARIO_ID==2038 || SCENARIO_ID==2048
            || SCENARIO_ID==4024){
            CO2S = 3; /* RCP6.0 */
            /**/
            if((fp_s[23]=fopen("./data/image_a1b_fcrop.dat","rt"))==NULL){
                printf("NO image_a1b_fcrop.dat !!\n");
                exit(1);
            }else{
                Flag_FOPEN[23] ++;
            }
            fp_s[45]=fopen("./data/image_a1b_fgrass.dat","rt");
        }
    }
    
    /* GEO-MIP */
    if(SCENARIO_ID >= 3000 && SCENARIO_ID <= 3999){
        CO2S = 2; /* RCP4.5 */
        /**/
        if((fp_s[23]=fopen("./data/image_a1b_fcrop.dat","rt"))==NULL){
			printf("NO image_a1b_fcrop.dat !!\n");
			exit(1);
		}else{
            Flag_FOPEN[23] ++;
        }
		fp_s[45]=fopen("./data/image_a1b_fgrass.dat","rt");
    }
    
    /* ISI-MIP 2 (historical) */
    /* if(SCENARIO_ID>=5000 && SCENARIO_ID<=5999){
        CO2S = 2;
        if((fp_s[23]=fopen("./data/image_a1b_fcrop.dat","rt"))==NULL){
			printf("NO image_a1b_fcrop.dat !!\n");
			exit(1);
		}
		fp_s[45]=fopen("./data/image_a1b_fgrass.dat","rt");
    } */
	
    /* ICARUS 2016/08/12 by A.Ito ***********************************/
    if(SCENARIO_ID==2201 || SCENARIO_ID==2203 || SCENARIO_ID==2208 || SCENARIO_ID==2210 ||
        SCENARIO_ID==2221 || SCENARIO_ID==2223 || SCENARIO_ID==2228 || SCENARIO_ID==2230 ||
        SCENARIO_ID==2241 || SCENARIO_ID==2243 || SCENARIO_ID==2248 || SCENARIO_ID==2250 ||
        SCENARIO_ID==2261 || SCENARIO_ID==2263 || SCENARIO_ID==2268 || SCENARIO_ID==2270 ||
        SCENARIO_ID==2281 || SCENARIO_ID==2283 || SCENARIO_ID==2288 || SCENARIO_ID==2290){
        CO2S = 1; /* RCP2.6 */
        /**/
        
        switch(SCENARIO_ID){
            case 2201:
                if((fp_s[23]=fopen("./data/fcrop_gfdl_aim_ssp1_rcp26_v1.txt","rt"))==NULL){
                    printf("NO fcrop_gfdl_aim_ssp1_rcp26_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2203:
                if((fp_s[23]=fopen("./data/fcrop_gfdl_aim_ssp2_rcp26_v1.txt","rt"))==NULL){
                    printf("NO fcrop_gfdl_aim_ssp2_rcp26_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2208:
                if((fp_s[23]=fopen("./data/fcrop_gfdl_aim_ssp4_rcp26_v1.txt","rt"))==NULL){
                    printf("NO fcrop_gfdl_aim_ssp4_rcp26_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2210:
                if((fp_s[23]=fopen("./data/fcrop_gfdl_aim_ssp5_rcp26_v1.txt","rt"))==NULL){
                    printf("NO fcrop_gfdl_aim_ssp5_rcp26_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2221:
                if((fp_s[23]=fopen("./data/fcrop_hadg_aim_ssp1_rcp26_v1.txt","rt"))==NULL){
                    printf("NO fcrop_hadg_aim_ssp1_rcp26_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2223:
                if((fp_s[23]=fopen("./data/fcrop_hadg_aim_ssp2_rcp26_v1.txt","rt"))==NULL){
                    printf("NO fcrop_hadg_aim_ssp2_rcp26_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2228:
                if((fp_s[23]=fopen("./data/fcrop_hadg_aim_ssp4_rcp26_v1.txt","rt"))==NULL){
                    printf("NO fcrop_hadg_aim_ssp4_rcp26_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2230:
                if((fp_s[23]=fopen("./data/fcrop_hadg_aim_ssp5_rcp26_v1.txt","rt"))==NULL){
                    printf("NO fcrop_hadg_aim_ssp5_rcp26_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2241:
                if((fp_s[23]=fopen("./data/fcrop_ipsl_aim_ssp1_rcp26_v1.txt","rt"))==NULL){
                    printf("NO fcrop_ipsl_aim_ssp1_rcp26_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2243:
                if((fp_s[23]=fopen("./data/fcrop_ipsl_aim_ssp2_rcp26_v1.txt","rt"))==NULL){
                    printf("NO fcrop_ipsl_aim_ssp2_rcp26_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2248:
                if((fp_s[23]=fopen("./data/fcrop_ipsl_aim_ssp4_rcp26_v1.txt","rt"))==NULL){
                    printf("NO fcrop_ipsl_aim_ssp4_rcp26_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2250:
                if((fp_s[23]=fopen("./data/fcrop_ipsl_aim_ssp5_rcp26_v1.txt","rt"))==NULL){
                    printf("NO fcrop_ipsl_aim_ssp5_rcp26_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2261:
                if((fp_s[23]=fopen("./data/fcrop_miro_aim_ssp1_rcp26_v1.txt","rt"))==NULL){
                    printf("NO fcrop_miro_aim_ssp1_rcp26_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2263:
                if((fp_s[23]=fopen("./data/fcrop_miro_aim_ssp2_rcp26_v1.txt","rt"))==NULL){
                    printf("NO fcrop_miro_aim_ssp2_rcp26_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2268:
                if((fp_s[23]=fopen("./data/fcrop_miro_aim_ssp4_rcp26_v1.txt","rt"))==NULL){
                    printf("NO fcrop_miro_aim_ssp4_rcp26_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2270:
                if((fp_s[23]=fopen("./data/fcrop_miro_aim_ssp5_rcp26_v1.txt","rt"))==NULL){
                    printf("NO fcrop_miro_aim_ssp5_rcp26_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2281:
                if((fp_s[23]=fopen("./data/fcrop_nore_aim_ssp1_rcp26_v1.txt","rt"))==NULL){
                    printf("NO fcrop_nore_aim_ssp1_rcp26_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2283:
                if((fp_s[23]=fopen("./data/fcrop_nore_aim_ssp2_rcp26_v1.txt","rt"))==NULL){
                    printf("NO fcrop_nore_aim_ssp2_rcp26_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2288:
                if((fp_s[23]=fopen("./data/fcrop_nore_aim_ssp4_rcp26_v1.txt","rt"))==NULL){
                    printf("NO fcrop_nore_aim_ssp4_rcp26_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2290:
                if((fp_s[23]=fopen("./data/fcrop_nore_aim_ssp5_rcp26_v1.txt","rt"))==NULL){
                    printf("NO fcrop_nore_aim_ssp5_rcp26_v1.txt !!\n");
                    exit(1);
                }
                break;
            default:
                printf("NO fcrop data !!\n");
                exit(1);
        }
        
        /* assumption: SRES data were used, instead of RCP data */
		fp_s[45]=fopen("./data/image_b1_fgrass.dat","rt");
    }else if(SCENARIO_ID==2202 || SCENARIO_ID==2204 || SCENARIO_ID==2206 || SCENARIO_ID==2209 || SCENARIO_ID==2211 ||
        SCENARIO_ID==2222 || SCENARIO_ID==2224 || SCENARIO_ID==2226 || SCENARIO_ID==2229 || SCENARIO_ID==2231 ||
        SCENARIO_ID==2242 || SCENARIO_ID==2244 || SCENARIO_ID==2246 || SCENARIO_ID==2249 || SCENARIO_ID==2251 ||
        SCENARIO_ID==2262 || SCENARIO_ID==2264 || SCENARIO_ID==2266 || SCENARIO_ID==2269 || SCENARIO_ID==2271 ||
        SCENARIO_ID==2282 || SCENARIO_ID==2284 || SCENARIO_ID==2286 || SCENARIO_ID==2289 || SCENARIO_ID==2291 ){
        CO2S = 2; /* RCP4.5 */
        /**/
        
        switch(SCENARIO_ID){
            case 2202:
                if((fp_s[23]=fopen("./data/fcrop_gfdl_aim_ssp1_rcp45_v1.txt","rt"))==NULL){
                    printf("NO fcrop_gfdl_aim_ssp1_rcp45_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2204:
                if((fp_s[23]=fopen("./data/fcrop_gfdl_aim_ssp2_rcp45_v1.txt","rt"))==NULL){
                    printf("NO fcrop_gfdl_aim_ssp2_rcp45_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2206:
                if((fp_s[23]=fopen("./data/fcrop_gfdl_aim_ssp3_rcp45_v1.txt","rt"))==NULL){
                    printf("NO fcrop_gfdl_aim_ssp3_rcp45_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2209:
                if((fp_s[23]=fopen("./data/fcrop_gfdl_aim_ssp4_rcp45_v1.txt","rt"))==NULL){
                    printf("NO fcrop_gfdl_aim_ssp4_rcp45_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2211:
                if((fp_s[23]=fopen("./data/fcrop_gfdl_aim_ssp5_rcp45_v1.txt","rt"))==NULL){
                    printf("NO fcrop_gfdl_aim_ssp5_rcp45_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2222:
                if((fp_s[23]=fopen("./data/fcrop_hadg_aim_ssp1_rcp45_v1.txt","rt"))==NULL){
                    printf("NO fcrop_hadg_aim_ssp1_rcp45_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2224:
                if((fp_s[23]=fopen("./data/fcrop_hadg_aim_ssp2_rcp45_v1.txt","rt"))==NULL){
                    printf("NO fcrop_hadg_aim_ssp2_rcp45_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2226:
                if((fp_s[23]=fopen("./data/fcrop_hadg_aim_ssp3_rcp45_v1.txt","rt"))==NULL){
                    printf("NO fcrop_hadg_aim_ssp3_rcp45_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2229:
                if((fp_s[23]=fopen("./data/fcrop_hadg_aim_ssp4_rcp45_v1.txt","rt"))==NULL){
                    printf("NO fcrop_hadg_aim_ssp4_rcp45_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2231:
                if((fp_s[23]=fopen("./data/fcrop_hadg_aim_ssp5_rcp45_v1.txt","rt"))==NULL){
                    printf("NO fcrop_hadg_aim_ssp5_rcp45_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2242:
                if((fp_s[23]=fopen("./data/fcrop_ipsl_aim_ssp1_rcp45_v1.txt","rt"))==NULL){
                    printf("NO fcrop_ipsl_aim_ssp1_rcp45_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2244:
                if((fp_s[23]=fopen("./data/fcrop_ipsl_aim_ssp2_rcp45_v1.txt","rt"))==NULL){
                    printf("NO fcrop_ipsl_aim_ssp2_rcp45_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2246:
                if((fp_s[23]=fopen("./data/fcrop_ipsl_aim_ssp3_rcp45_v1.txt","rt"))==NULL){
                    printf("NO fcrop_ipsl_aim_ssp3_rcp45_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2249:
                if((fp_s[23]=fopen("./data/fcrop_ipsl_aim_ssp4_rcp45_v1.txt","rt"))==NULL){
                    printf("NO fcrop_ipsl_aim_ssp4_rcp45_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2251:
                if((fp_s[23]=fopen("./data/fcrop_ipsl_aim_ssp5_rcp45_v1.txt","rt"))==NULL){
                    printf("NO fcrop_ipsl_aim_ssp5_rcp45_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2262:
                if((fp_s[23]=fopen("./data/fcrop_miro_aim_ssp1_rcp45_v1.txt","rt"))==NULL){
                    printf("NO fcrop_miro_aim_ssp1_rcp45_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2264:
                if((fp_s[23]=fopen("./data/fcrop_miro_aim_ssp2_rcp45_v1.txt","rt"))==NULL){
                    printf("NO fcrop_miro_aim_ssp2_rcp45_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2266:
                if((fp_s[23]=fopen("./data/fcrop_miro_aim_ssp3_rcp45_v1.txt","rt"))==NULL){
                    printf("NO fcrop_miro_aim_ssp3_rcp45_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2269:
                if((fp_s[23]=fopen("./data/fcrop_miro_aim_ssp4_rcp45_v1.txt","rt"))==NULL){
                    printf("NO fcrop_miro_aim_ssp4_rcp45_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2271:
                if((fp_s[23]=fopen("./data/fcrop_miro_aim_ssp5_rcp45_v1.txt","rt"))==NULL){
                    printf("NO fcrop_miro_aim_ssp5_rcp45_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2282:
                if((fp_s[23]=fopen("./data/fcrop_nore_aim_ssp1_rcp45_v1.txt","rt"))==NULL){
                    printf("NO fcrop_nore_aim_ssp1_rcp45_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2284:
                if((fp_s[23]=fopen("./data/fcrop_nore_aim_ssp2_rcp45_v1.txt","rt"))==NULL){
                    printf("NO fcrop_nore_aim_ssp2_rcp45_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2286:
                if((fp_s[23]=fopen("./data/fcrop_nore_aim_ssp3_rcp45_v1.txt","rt"))==NULL){
                    printf("NO fcrop_nore_aim_ssp3_rcp45_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2289:
                if((fp_s[23]=fopen("./data/fcrop_nore_aim_ssp4_rcp45_v1.txt","rt"))==NULL){
                    printf("NO fcrop_nore_aim_ssp4_rcp45_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2291:
                if((fp_s[23]=fopen("./data/fcrop_nore_aim_ssp5_rcp45_v1.txt","rt"))==NULL){
                    printf("NO fcrop_nore_aim_ssp5_rcp45_v1.txt !!\n");
                    exit(1);
                }
                break;
            default:
                printf("NO fcrop data !!\n");
                exit(1);
        }
        
        /* assumption: SRES data were used, instead of RCP data */
		fp_s[45]=fopen("./data/image_a1b_fgrass.dat","rt");
    }else if(SCENARIO_ID == 2205 || SCENARIO_ID == 2207 ||
        SCENARIO_ID==2225 || SCENARIO_ID==2227 ||
        SCENARIO_ID==2245 || SCENARIO_ID==2247 ||
        SCENARIO_ID==2265 || SCENARIO_ID==2267 ||
        SCENARIO_ID==2285 || SCENARIO_ID==2287){
        CO2S = 3; /* RCP6.0 */
        /**/
        switch(SCENARIO_ID){
            case 2205:
                if((fp_s[23]=fopen("./data/fcrop_gfdl_aim_ssp2_rcp60_v1.txt","rt"))==NULL){
                    printf("NO fcrop_gfdl_aim_ssp2_rcp60_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2207:
                if((fp_s[23]=fopen("./data/fcrop_gfdl_aim_ssp3_rcp60_v1.txt","rt"))==NULL){
                    printf("NO fcrop_gfdl_aim_ssp3_rcp60_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2225:
                if((fp_s[23]=fopen("./data/fcrop_hadg_aim_ssp2_rcp60_v1.txt","rt"))==NULL){
                    printf("NO fcrop_hadg_aim_ssp2_rcp60_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2227:
                if((fp_s[23]=fopen("./data/fcrop_hadg_aim_ssp3_rcp60_v1.txt","rt"))==NULL){
                    printf("NO fcrop_hadg_aim_ssp3_rcp60_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2245:
                if((fp_s[23]=fopen("./data/fcrop_ipsl_aim_ssp2_rcp60_v1.txt","rt"))==NULL){
                    printf("NO fcrop_ipsl_aim_ssp2_rcp60_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2247:
                if((fp_s[23]=fopen("./data/fcrop_ipsl_aim_ssp3_rcp60_v1.txt","rt"))==NULL){
                    printf("NO fcrop_ipsl_aim_ssp3_rcp60_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2265:
                if((fp_s[23]=fopen("./data/fcrop_miro_aim_ssp2_rcp60_v1.txt","rt"))==NULL){
                    printf("NO fcrop_miro_aim_ssp2_rcp60_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2267:
                if((fp_s[23]=fopen("./data/fcrop_miro_aim_ssp3_rcp60_v1.txt","rt"))==NULL){
                    printf("NO fcrop_miro_aim_ssp3_rcp60_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2285:
                if((fp_s[23]=fopen("./data/fcrop_nore_aim_ssp2_rcp60_v1.txt","rt"))==NULL){
                    printf("NO fcrop_nore_aim_ssp2_rcp60_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2287:
                if((fp_s[23]=fopen("./data/fcrop_nore_aim_ssp3_rcp60_v1.txt","rt"))==NULL){
                    printf("NO fcrop_nore_aim_ssp3_rcp60_v1.txt !!\n");
                    exit(1);
                }
                break;
            default:
                printf("NO fcrop data !!\n");
                exit(1);
        }
        
        /* assumption: SRES data were used, instead of RCP data */
		if((fp_s[45]=fopen("./data/image_a1b_fgrass.dat","rt"))==NULL){
 			printf("NO image_a1b_fgrass.dat !!\n");
			exit(1);
        }else{
            Flag_FOPEN[45] = 1;
        }
    }else if(SCENARIO_ID == 2601){
        CO2S = 0; /* fixed */
        if((fp_s[23]=fopen("./data/fcrop_ssp1_current_v1.txt","rt"))==NULL){
            printf("NO fcrop_ssp1_current_v1.txt !!\n");
            exit(1);
        }else{
            Flag_FOPEN[23] ++;
        }
        /* assumption: SRES data were used, instead of RCP data */
		if((fp_s[45]=fopen("./data/image_a1b_fgrass.dat","rt"))==NULL){
 			printf("NO image_a1b_fgrass.dat !!\n");
			exit(1);
        }else{
            Flag_FOPEN[45] = 1;
        }
    }else if(SCENARIO_ID == 2602){
        CO2S = 0; /* fixed */
        if((fp_s[23]=fopen("./data/fcrop_ssp2_current_v1.txt","rt"))==NULL){
            printf("NO fcrop_ssp2_current_v1.txt !!\n");
            exit(1);
        }else{
            Flag_FOPEN[23] ++;
        }
        /* assumption: SRES data were used, instead of RCP data */
		if((fp_s[45]=fopen("./data/image_a1b_fgrass.dat","rt"))==NULL){
 			printf("NO image_a1b_fgrass.dat !!\n");
			exit(1);
        }else{
            Flag_FOPEN[45] = 1;
        }
    }else if(SCENARIO_ID == 2603){
        CO2S = 0; /* fixed */
        if((fp_s[23]=fopen("./data/fcrop_ssp3_current_v1.txt","rt"))==NULL){
            printf("NO fcrop_ssp3_current_v1.txt !!\n");
            exit(1);
        }
        /* assumption: SRES data were used, instead of RCP data */
		if((fp_s[45]=fopen("./data/image_a1b_fgrass.dat","rt"))==NULL){
 			printf("NO image_a1b_fgrass.dat !!\n");
			exit(1);
        }else{
            Flag_FOPEN[45] = 1;
        }
    }
    
    /* GEO-MIP ***************************************/
    if(SCENARIO_ID>=3000 && SCENARIO_ID<=3999){
        CO2S = 2; /* RCP4.5 */
        /**/
        if((fp_s[23]=fopen("./data/image_a1b_fcrop.dat","rt"))==NULL){
			printf("NO image_a1b_fcrop.dat !!\n");
			exit(1);
		}else{
            Flag_FOPEN[23] ++;
        }
		if((fp_s[45]=fopen("./data/image_a1b_fgrass.dat","rt"))==NULL){
 			printf("NO image_a1b_fgrass.dat !!\n");
			exit(1);
        }else{
            Flag_FOPEN[45] ++;
        }
    }
    
    /* ISI-MIP2 ******************************/
    if(SCENARIO_ID>=5000 && SCENARIO_ID<=5009){
        /* ISI-MIP2a (historical) */
        CO2S = 2; /* RCP4.5 */
        /**/
        if((fp_s[23]=fopen("./data/image_a1b_fcrop.dat","rt"))==NULL){
			printf("NO image_a1b_fcrop.dat !!\n");
			exit(1);
		}else{
            Flag_FOPEN[23] ++;
        }
		if((fp_s[45]=fopen("./data/image_a1b_fgrass.dat","rt"))==NULL){
 			printf("NO image_a1b_fgrass.dat !!\n");
			exit(1);
        }else{
            Flag_FOPEN[45] ++;
        }
    }else if(SCENARIO_ID>=5010 && SCENARIO_ID<=5100){
        /* ISI-MIP2b (1.5/2.0deg) */
        
        if(SCENARIO_ID == 5010 || SCENARIO_ID == 5020 || SCENARIO_ID == 5030 || SCENARIO_ID == 5040){
            CO2S = 1;
        }
        if(SCENARIO_ID == 5011 || SCENARIO_ID == 5021 || SCENARIO_ID == 5031 || SCENARIO_ID == 5041){
            CO2S = 2;
        }
        if(SCENARIO_ID == 5012 || SCENARIO_ID == 5022 || SCENARIO_ID == 5032 || SCENARIO_ID == 5042){
            CO2S = 3;
        }
        if(SCENARIO_ID == 5013 || SCENARIO_ID == 5023 || SCENARIO_ID == 5033 || SCENARIO_ID == 5043){
            CO2S = 1;
        }
        if(SCENARIO_ID == 5014 || SCENARIO_ID == 5024 || SCENARIO_ID == 5034 || SCENARIO_ID == 5044){
            CO2S = 1;
        }
        if(SCENARIO_ID == 5015 || SCENARIO_ID == 5025 || SCENARIO_ID == 5035 || SCENARIO_ID == 5045){
            CO2S = 2;
        }
        if(SCENARIO_ID == 5016 || SCENARIO_ID == 5026 || SCENARIO_ID == 5036 || SCENARIO_ID == 5046){
            CO2S = 3;
        }
        if(SCENARIO_ID == 5060 || SCENARIO_ID == 5070 || SCENARIO_ID == 5080 || SCENARIO_ID == 5090){
            CO2S = 1;
        }
        if(SCENARIO_ID == 5061 || SCENARIO_ID == 5071 || SCENARIO_ID == 5081 || SCENARIO_ID == 5091){
            CO2S = 2;
        }
        if(SCENARIO_ID == 5062 || SCENARIO_ID == 5072 || SCENARIO_ID == 5082 || SCENARIO_ID == 5092){
            CO2S = 3;
        }
        if(SCENARIO_ID == 5063 || SCENARIO_ID == 5073 || SCENARIO_ID == 5083 || SCENARIO_ID == 5093){
            CO2S = 3;
        }
        if(SCENARIO_ID == 5017 || SCENARIO_ID == 5027 || SCENARIO_ID == 5037 || SCENARIO_ID == 5047){
            /* RCP8.5: added 2019/02/25 by A.Ito */
            CO2S = 4;
        }

        /* dummy: won't be used */
        if((fp_s[23]=fopen("./data/image_a1b_fcrop.dat","rt"))==NULL){
            printf("NO image_a1b_fcrop.dat !!\n");
            exit(1);
        }else{
            Flag_FOPEN[23] ++;
        }
        if((fp_s[45]=fopen("./data/image_a1b_fgrass.dat","rt"))==NULL){
             printf("NO image_a1b_fgrass.dat !!\n");
            exit(1);
        }else{
            Flag_FOPEN[45] ++;
        }
    }

   /**************************************************************/
    /* SAGE crop data */
	if( (fp_s[46]=fopen("./data/sage_crops.dat","rt"))==NULL ){  
		printf("No sage_crops.dat\n");  
		exit(1); 
	}else{
        Flag_FOPEN[46] ++;
    }
    
    /* radiation diffusion parameters */
	if( (fp_s[47]=fopen("./data/reg_raddif_hd.dat","rt"))==NULL ){  
		printf("No diffuse rad model data\n");  
		exit(1); 
	}else{
        Flag_FOPEN[47] ++;
    }
	
    /* SSM/I inundation data */
	if( (fp_s[48]=fopen("./data/ssmi_season.dat","rt"))==NULL ){  
		printf("No SSMI inundation data\n");  
		exit(1); 
	}else{
        Flag_FOPEN[48] ++;
    }
    
    /* permafrost data */
	if( (fp_s[49]=fopen("./data/permafrost_nsidc.dat","rt"))==NULL ){  
		printf("No permafrost data\n");  
		exit(1); 
	}else{
        Flag_FOPEN[49] ++;
    }
	
	/* crop and pasture fractions: 1700-2007 */
	/* Ramankutty & Kimball: added 2010/07/07 (A.Ito) */
	if( (fp_s[50]=fopen("./data/glcrop_1700-2007_0.5.dat","rt"))==NULL ){  
		printf("No glcrop_1700-2007_0.5.dat\n");  
		exit(1); 
	}else{
        Flag_FOPEN[50] ++;
    }
	
	if( (fp_s[51]=fopen("./data/glpast_1700-2007_0.5.dat","rt"))==NULL ){  
		printf("No glpast_1700-2007_0.5.dat\n");  
		exit(1); 
	}else{
        Flag_FOPEN[51] ++;
    }
	
    /* N deposition **************************************************/
    if(SCENARIO_ID==5010 || SCENARIO_ID==5020 || SCENARIO_ID==5030 || SCENARIO_ID==5040 ||
        SCENARIO_ID==5060 || SCENARIO_ID==5070 || SCENARIO_ID==5080 || SCENARIO_ID==5090){
        if( (fp_s[25]=fopen("./data/ndepo_ann_2005soc_1661-2299.flt","rb"))==NULL ){
            printf("No ndepo_ann_2005soc_1661-2299.flt\n");
            exit(1); 
        }else{
            Flag_FOPEN[25] ++;
        }
    }else if(SCENARIO_ID==5011 || SCENARIO_ID==5021 || SCENARIO_ID==5031 || SCENARIO_ID==5041 ||
            SCENARIO_ID==5013 || SCENARIO_ID==5023 || SCENARIO_ID==5033 || SCENARIO_ID==5043 ||
            SCENARIO_ID==5015 || SCENARIO_ID==5025 || SCENARIO_ID==5035 || SCENARIO_ID==5045 ||
            SCENARIO_ID==5061 || SCENARIO_ID==5071 || SCENARIO_ID==5081 || SCENARIO_ID==5091){
        if( (fp_s[25]=fopen("./data/ndepo_ann_rcp26_1661-2299.flt","rb"))==NULL ){
            printf("No ndepo_ann_rcp26_1661-2299.flt\n");
            exit(1); 
        }else{
            Flag_FOPEN[25] ++;
        }
    }else if(SCENARIO_ID==5012 || SCENARIO_ID==5022 || SCENARIO_ID==5032 || SCENARIO_ID==5042 ||
            SCENARIO_ID==5014 || SCENARIO_ID==5024 || SCENARIO_ID==5034 || SCENARIO_ID==5044 ||
            SCENARIO_ID==5016 || SCENARIO_ID==5026 || SCENARIO_ID==5036 || SCENARIO_ID==5046 ||
            SCENARIO_ID==5062 || SCENARIO_ID==5072 || SCENARIO_ID==5082 || SCENARIO_ID==5092 ||
            SCENARIO_ID==5063 || SCENARIO_ID==5073 || SCENARIO_ID==5083 || SCENARIO_ID==5093){
        if( (fp_s[25]=fopen("./data/ndepo_ann_rcp60_1661-2299.flt","rb"))==NULL ){
            printf("No ndepo_ann_rcp60_1661-2299.flt\n");
            exit(1); 
        }else{
            Flag_FOPEN[25] ++;
        }
    }else if(SCENARIO_ID==5017 || SCENARIO_ID==5027 || SCENARIO_ID==5037 || SCENARIO_ID==5047){
        if( (fp_s[25]=fopen("./data/ndepo_ann_rcp85_1661-2299.flt","rb"))==NULL ){
            printf("No ndepo_ann_rcp85_1661-2299.flt\n");
            exit(1);
        }else{
            Flag_FOPEN[25] ++;
        }
    }else{
        /* by Galloway & Dentener */
        if( (fp_s[25]=fopen("./data/N-deposition.dat","rt"))==NULL ){
            printf("No N-deposition.dat\n");  
            exit(1); 
        }else{
            Flag_FOPEN[25] ++;
        }
    }
    
	/* by CHASER */
	if( (fp_s[52]=fopen("./data/chaser_ndepo_2001mon.txt","rt"))==NULL ){  
		printf("No CHASER N deposition data: chaser_ndepo_2001mon.txt\n");
		exit(1); 
	}else{
        Flag_FOPEN[52] ++;
    }
	if( (fp_s[85]=fopen("./data/CHASER_V4.0+HTAP2_preind-1851-1855.txt","rt"))==NULL ){
		printf("No CHASER N deposition data: CHASER_V4.0+HTAP2_preind-1851-1855.txt\n");
		exit(1); 
	}else{
        Flag_FOPEN[85] ++;
    }
	if( (fp_s[86]=fopen("./data/CHASER_V4.0+CCMI_REFC1SD_2006-2010.txt","rt"))==NULL ){
		printf("No CHASER N deposition data: CHASER_V4.0+CCMI_REFC1SD_2006-2010.txt\n");
		exit(1); 
	}else{
        Flag_FOPEN[86] ++;
    }
	
	/* wood harvest based on RCP-harmonized data: LUHa.v1 *************/
	/* added by A.Ito (2010/10/15) */
    if(LANDUSE == 9 || LANDUSE == 10 || LANDUSE == 11 || LANDUSE == 12 ||
        LANDUSE == 13 || LANDUSE == 17){
        if( (fp_s[53]=fopen("./data/luh_gvbh1_1500_2005.txt","rt"))==NULL ){  
            printf("No luh_gvbh1_1500_2005.txt\n");  
            exit(1); 
        }else{
            Flag_FOPEN[53] ++;
        }
        if( (fp_s[54]=fopen("./data/luh_gvbh2_1500_2005.txt","rt"))==NULL ){  
            printf("No luh_gvbh2_1500_2005.txt\n");  
            exit(1); 
        }else{
            Flag_FOPEN[54] ++;
        }
        if( (fp_s[55]=fopen("./data/luh_gsbh1_1500_2005.txt","rt"))==NULL ){  
            printf("No luh_gsbh1_1500_2005.txt\n");  
            exit(1); 
        }else{
            Flag_FOPEN[55] ++;
        }
        if( (fp_s[56]=fopen("./data/luh_gsbh2_1500_2005.txt","rt"))==NULL ){  
            printf("No luh_gsbh2_1500_2005.txt\n");  
            exit(1); 
        }else{
            Flag_FOPEN[56] ++;
        }
        if( (fp_s[57]=fopen("./data/luh_gsbh3_1500_2005.txt","rt"))==NULL ){  
            printf("No luh_gsbh3_1500_2005.txt\n");  
            exit(1); 
        }else{
            Flag_FOPEN[57] ++;
        }
    }else{
        if( (fp_s[53]=fopen("./data/luc_eos2_gvbh1_1700-2005.dat","rt"))==NULL ){  
            printf("No luc_eos2_gvbh1_1700-2005.dat\n");  
            exit(1); 
        }
        if( (fp_s[54]=fopen("./data/luc_eos2_gvbh2_1700-2005.dat","rt"))==NULL ){  
            printf("No luc_eos2_gvbh2_1700-2005.dat\n");  
            exit(1); 
        }
        if( (fp_s[55]=fopen("./data/luc_eos2_hsbh1_1700-2005.dat","rt"))==NULL ){  
            printf("No luc_eos2_hsbh1_1700-2005.dat\n");  
            exit(1); 
        }
        if( (fp_s[56]=fopen("./data/luc_eos2_hsbh2_1700-2005.dat","rt"))==NULL ){  
            printf("No luc_eos2_hsbh2_1700-2005.dat\n");  
            exit(1); 
        }
        if( (fp_s[57]=fopen("./data/luc_eos2_hsbh3_1700-2005.dat","rt"))==NULL ){  
            printf("No luc_eos2_hsbh3_1700-2005.dat\n");  
            exit(1); 
        }else{
            Flag_FOPEN[57] ++;
        }
    }
        
    if(LANDUSE == 9 || LANDUSE == 10){
        if( (fp_s[78]=fopen("./data/luh_gvbh1_rcp45_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gvbh1_rcp45_2005_2100.txt\n");  
            exit(1); 
        }
        if( (fp_s[79]=fopen("./data/luh_gvbh2_rcp45_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gvbh2_rcp45_2005_2100.txt\n");  
            exit(1); 
        }
        if( (fp_s[80]=fopen("./data/luh_gsbh1_rcp45_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gsbh1_rcp45_2005_2100.txt\n");  
            exit(1); 
        }
        if( (fp_s[81]=fopen("./data/luh_gsbh2_rcp45_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gsbh2_rcp45_2005_2100.txt\n");  
            exit(1); 
        }
        if( (fp_s[82]=fopen("./data/luh_gsbh3_rcp45_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gsbh3_rcp45_2005_2100.txt\n");  
            exit(1); 
        }
    }else if(LANDUSE == 11 || LANDUSE == 17){
        if( (fp_s[78]=fopen("./data/luh_gvbh1_rcp26_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gvbh1_rcp26_2005_2100.txt\n");  
            exit(1); 
        }
        if( (fp_s[79]=fopen("./data/luh_gvbh2_rcp26_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gvbh2_rcp26_2005_2100.txt\n");  
            exit(1); 
        }
        if( (fp_s[80]=fopen("./data/luh_gsbh1_rcp26_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gsbh1_rcp26_2005_2100.txt\n");  
            exit(1); 
        }
        if( (fp_s[81]=fopen("./data/luh_gsbh2_rcp26_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gsbh2_rcp26_2005_2100.txt\n");  
            exit(1); 
        }
        if( (fp_s[82]=fopen("./data/luh_gsbh3_rcp26_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gsbh3_rcp26_2005_2100.txt\n");  
            exit(1); 
        }
    }else if(LANDUSE == 12){
        if( (fp_s[78]=fopen("./data/luh_gvbh1_rcp60_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gvbh1_rcp60_2005_2100.txt\n");  
            exit(1); 
        }
        if( (fp_s[79]=fopen("./data/luh_gvbh2_rcp60_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gvbh2_rcp60_2005_2100.txt\n");  
            exit(1); 
        }
        if( (fp_s[80]=fopen("./data/luh_gsbh1_rcp60_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gsbh1_rcp60_2005_2100.txt\n");  
            exit(1); 
        }
        if( (fp_s[81]=fopen("./data/luh_gsbh2_rcp60_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gsbh2_rcp60_2005_2100.txt\n");  
            exit(1); 
        }
        if( (fp_s[82]=fopen("./data/luh_gsbh3_rcp60_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gsbh3_rcp60_2005_2100.txt\n");  
            exit(1); 
        }
    }else if(LANDUSE == 13){
        if( (fp_s[78]=fopen("./data/luh_gvbh1_rcp85_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gvbh1_rcp85_2005_2100.txt\n");  
            exit(1); 
        }
        if( (fp_s[79]=fopen("./data/luh_gvbh2_rcp85_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gvbh2_rcp85_2005_2100.txt\n");  
            exit(1); 
        }
        if( (fp_s[80]=fopen("./data/luh_gsbh1_rcp85_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gsbh1_rcp85_2005_2100.txt\n");  
            exit(1); 
        }
        if( (fp_s[81]=fopen("./data/luh_gsbh2_rcp85_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gsbh2_rcp85_2005_2100.txt\n");  
            exit(1); 
        }
        if( (fp_s[82]=fopen("./data/luh_gsbh3_rcp85_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gsbh3_rcp85_2005_2100.txt\n");  
            exit(1); 
        }
    }else if(LANDUSE == 14){
        /* ICARUS SSP data by T.Kinoshita: 2014/09/04 => 09/18 ******************/
        if( (fp_s[78]=fopen("./data/luc_wh_ssp1_2001-2100.txt","rt"))==NULL ){
            printf("No luc_wh_ssp1_2001-2100.txt\n");
            exit(1); 
        }
        if( (fp_s[79]=fopen("./data/luh_gvbh2_rcp45_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gvbh2_rcp45_2005_2100.txt\n");  
            exit(1); 
        }
        if( (fp_s[80]=fopen("./data/luh_gsbh1_rcp45_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gsbh1_rcp45_2005_2100.txt\n");  
            exit(1); 
        }
        if( (fp_s[81]=fopen("./data/luh_gsbh2_rcp45_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gsbh2_rcp45_2005_2100.txt\n");  
            exit(1); 
        }
        if( (fp_s[82]=fopen("./data/luh_gsbh3_rcp45_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gsbh3_rcp45_2005_2100.txt\n");  
            exit(1); 
        }
    }else if(LANDUSE == 15){
        /* ICARUS SSP data by T.Kinoshita: 2014/09/04 => 09/18 ******************/
        if( (fp_s[78]=fopen("./data/luc_wh_ssp2_2001-2100.txt","rt"))==NULL ){
            printf("No luc_wh_ssp2_2001-2100.txt\n");
            exit(1); 
        }else{
            Flag_FOPEN[78] ++;
        }
        if( (fp_s[79]=fopen("./data/luh_gvbh2_rcp45_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gvbh2_rcp45_2005_2100.txt\n");  
            exit(1); 
        }else{
            Flag_FOPEN[79] ++;
        }
        if( (fp_s[80]=fopen("./data/luh_gsbh1_rcp45_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gsbh1_rcp45_2005_2100.txt\n");  
            exit(1); 
        }else{
            Flag_FOPEN[80] ++;
        }
        if( (fp_s[81]=fopen("./data/luh_gsbh2_rcp45_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gsbh2_rcp45_2005_2100.txt\n");  
            exit(1); 
        }else{
            Flag_FOPEN[81] ++;
        }
        if( (fp_s[82]=fopen("./data/luh_gsbh3_rcp45_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gsbh3_rcp45_2005_2100.txt\n");  
            exit(1); 
        }else{
            Flag_FOPEN[82] ++;
        }
    }else if(LANDUSE == 16){
        /* ICARUS SSP data by T.Kinoshita: 2014/09/04 => 09/18 ******************/
        if( (fp_s[78]=fopen("./data/luc_wh_ssp3_2001-2100.txt","rt"))==NULL ){
            printf("No luc_wh_ssp3_2001-2100.txt\n");
            exit(1); 
        }else{
            Flag_FOPEN[78] ++;
        }
        if( (fp_s[79]=fopen("./data/luh_gvbh2_rcp45_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gvbh2_rcp45_2005_2100.txt\n");  
            exit(1); 
        }else{
            Flag_FOPEN[79] ++;
        }
        if( (fp_s[80]=fopen("./data/luh_gsbh1_rcp45_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gsbh1_rcp45_2005_2100.txt\n");  
            exit(1); 
        }else{
            Flag_FOPEN[80] ++;
        }
        if( (fp_s[81]=fopen("./data/luh_gsbh2_rcp45_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gsbh2_rcp45_2005_2100.txt\n");  
            exit(1); 
        }else{
            Flag_FOPEN[81] ++;
        }
        if( (fp_s[82]=fopen("./data/luh_gsbh3_rcp45_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gsbh3_rcp45_2005_2100.txt\n");  
            exit(1); 
        }else{
            Flag_FOPEN[82] ++;
        }
    }else{
        if( (fp_s[78]=fopen("./data/luh_gvbh1_rcp45_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gvbh1_rcp45_2005_2100.txt\n");  
            exit(1); 
        }else{
            Flag_FOPEN[78] ++;
        }
        if( (fp_s[79]=fopen("./data/luh_gvbh2_rcp45_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gvbh2_rcp45_2005_2100.txt\n");  
            exit(1); 
        }else{
            Flag_FOPEN[79] ++;
        }
        if( (fp_s[80]=fopen("./data/luh_gsbh1_rcp45_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gsbh1_rcp45_2005_2100.txt\n");  
            exit(1); 
        }else{
            Flag_FOPEN[80] ++;
        }
        if( (fp_s[81]=fopen("./data/luh_gsbh2_rcp45_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gsbh2_rcp45_2005_2100.txt\n");  
            exit(1); 
        }else{
            Flag_FOPEN[81] ++;
        }
        if( (fp_s[82]=fopen("./data/luh_gsbh3_rcp45_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gsbh3_rcp45_2005_2100.txt\n");  
            exit(1); 
        }else{
            Flag_FOPEN[82] ++;
        }
    }
    
    if(ISIMIP_RUN == 4){
        /* ISI-MIP2b: 2016/12/24 by A.Ito */
        /* fp_s[87]=fopen("./data/landuse_ann_soc_1661-2299.flt","rb") */
        /* ISI-MIP2b: 2017/11/01, 2018/03/01 by A.Ito */

        if(LANDUSE == 25 || SCENARIO_ID==5010 || SCENARIO_ID==5020
                   || SCENARIO_ID==5030 || SCENARIO_ID==5040){
            if( (fp_s[87]=fopen("./data/landuse_ann_1860_hist_2005_2005_1661-2299.flt","rb"))==NULL ){
                printf("No landuse_ann_1860_hist_2005_2005_1661-2299.flt\n");
                exit(1);
            }else{
                Flag_FOPEN[87] ++;
            }
        }else if(LANDUSE == 24){
            if(SCENARIO_ID==5060 || SCENARIO_ID==5061){
                if( (fp_s[87]=fopen("./data/landuse_ann_1860_1860_1860_1860_1661-2299.flt","rb"))==NULL ){
                    printf("No landuse_ann_1860_1860_1860_1860_1661-2299.flt\n");
                    exit(1);
                }else{
                    Flag_FOPEN[87] ++;
                }
            }else if(SCENARIO_ID==5011 || SCENARIO_ID==5017){
                if( (fp_s[87]=fopen("./data/landuse_ann_1860_hist_2005_2005_1661-2299.flt","rb"))==NULL ){
                    printf("No landuse_ann_1860_hist_2005_2005_1661-2299.flt\n");
                    exit(1);
                }else{
                    Flag_FOPEN[87] ++;
                }
            }else if(SCENARIO_ID==5013 || SCENARIO_ID==5015){
                if( (fp_s[87]=fopen("./data/landuse_ann_1860_hist_rcp26_2100rcp26_1661-2299_ipsl.flt","rb"))==NULL ){
                    printf("No landuse_ann_1860_hist_rcp26_2100rcp26_1661-2299_ipsl.flt\n");
                    exit(1);
                }else{
                    Flag_FOPEN[87] ++;
                }
            }else if(SCENARIO_ID==5012 || SCENARIO_ID==5062){
                if( (fp_s[87]=fopen("./data/landuse_ann_1860_hist_2005_2005_1661-2299.flt","rb"))==NULL ){
                    printf("No landuse_ann_1860_hist_2005_2005_1661-2299.flt\n");
                    exit(1);
                }else{
                    Flag_FOPEN[87] ++;
                }
            }else if(SCENARIO_ID==5014 || SCENARIO_ID==5016){
                if( (fp_s[87]=fopen("./data/landuse_ann_1860_hist_rcp60_2100_1661-2299_ipsl.flt","rb"))==NULL ){
                    printf("No landuse_ann_1860_hist_rcp60_2100_1661-2299_ipsl.flt\n");
                    exit(1);
                }else{
                    Flag_FOPEN[87] ++;
                }
            }else if(SCENARIO_ID==5063){
                if( (fp_s[87]=fopen("./data/landuse_ann_1860_hist_1860_1860_1661-2299.flt","rb"))==NULL ){
                    printf("No landuse_ann_1860_hist_1860_1860_1661-2299.flt\n");
                    exit(1);
                }else{
                    Flag_FOPEN[87] ++;
                }
            }
            
            if(SCENARIO_ID==5070 || SCENARIO_ID==5071){
                if( (fp_s[87]=fopen("./data/landuse_ann_1860_1860_1860_1860_1661-2299.flt","rb"))==NULL ){
                    printf("No landuse_ann_1860_1860_1860_1860_1661-2299.flt\n");
                    exit(1);
                }else{
                    Flag_FOPEN[87] ++;
                }
            }else if(SCENARIO_ID==5021 || SCENARIO_ID==5027){
                if( (fp_s[87]=fopen("./data/landuse_ann_1860_hist_2005_2005_1661-2299.flt","rb"))==NULL ){
                    printf("No landuse_ann_1860_hist_2005_2005_1661-2299.flt\n");
                    exit(1);
                }else{
                    Flag_FOPEN[87] ++;
                }
            }else if(SCENARIO_ID==5023 || SCENARIO_ID==5025){
                if( (fp_s[87]=fopen("./data/landuse_ann_1860_hist_rcp26_2100rcp26_1661-2299_gfdl.flt","rb"))==NULL ){
                    printf("No landuse_ann_1860_hist_rcp26_2100rcp26_1661-2299_gfdl.flt\n");
                    exit(1);
                }else{
                    Flag_FOPEN[87] ++;
                }
            }else if(SCENARIO_ID==5022 || SCENARIO_ID==5072){
                if( (fp_s[87]=fopen("./data/landuse_ann_1860_hist_2005_2005_1661-2299.flt","rb"))==NULL ){
                    printf("No landuse_ann_1860_hist_2005_2005_1661-2299.flt\n");
                    exit(1);
                }else{
                    Flag_FOPEN[87] ++;
                }
            }else if(SCENARIO_ID==5024 || SCENARIO_ID==5026){
                if( (fp_s[87]=fopen("./data/landuse_ann_1860_hist_rcp60_2100_1661-2299_gfdl.flt","rb"))==NULL ){
                    printf("No landuse_ann_1860_hist_2005_2005_rcp60_2100_gfdl.flt\n");
                    exit(1);
                }else{
                    Flag_FOPEN[87] ++;
                }
            }else if(SCENARIO_ID==5073){
                if( (fp_s[87]=fopen("./data/landuse_ann_1860_hist_1860_1860_1661-2299.flt","rb"))==NULL ){
                    printf("No landuse_ann_1860_hist_1860_1860_1661-2299.flt\n");
                    exit(1);
                }else{
                    Flag_FOPEN[87] ++;
                }
            }

            if(SCENARIO_ID==5080 || SCENARIO_ID==5081){
                if( (fp_s[87]=fopen("./data/landuse_ann_1860_1860_1860_1860_1661-2299.flt","rb"))==NULL ){
                    printf("No landuse_ann_1860_1860_1860_1860_1661-2299.flt\n");
                    exit(1);
                }else{
                    Flag_FOPEN[87] ++;
                }
            }else if(SCENARIO_ID==5031 || SCENARIO_ID==5037){
                if( (fp_s[87]=fopen("./data/landuse_ann_1860_hist_2005_2005_1661-2299.flt","rb"))==NULL ){
                    printf("No landuse_ann_1860_hist_2005_2005_1661-2299.flt\n");
                    exit(1);
                }else{
                    Flag_FOPEN[87] ++;
                }
            }else if(SCENARIO_ID==5033 || SCENARIO_ID==5035){
                if( (fp_s[87]=fopen("./data/landuse_ann_1860_hist_rcp26_2100rcp26_1661-2299_miroc.flt","rb"))==NULL ){
                    printf("No landuse_ann_1860_hist_rcp26_2100rcp26_1661-2299_miroc.flt\n");
                    exit(1);
                }else{
                    Flag_FOPEN[87] ++;
                }
            }else if(SCENARIO_ID==5032 || SCENARIO_ID==5082){
                if( (fp_s[87]=fopen("./data/landuse_ann_1860_hist_2005_2005_1661-2299.flt","rb"))==NULL ){
                    printf("No landuse_ann_1860_hist_2005_2005_1661-2299.flt\n");
                    exit(1);
                }else{
                    Flag_FOPEN[87] ++;
                }
            }else if(SCENARIO_ID==5034 || SCENARIO_ID==5036){
                if( (fp_s[87]=fopen("./data/landuse_ann_1860_hist_rcp60_2100_1661-2299_miroc.flt","rb"))==NULL ){
                    printf("No landuse_ann_1860_hist_2005_2005_rcp60_2100_miroc.flt\n");
                    exit(1);
                }else{
                    Flag_FOPEN[87] ++;
                }
            }else if(SCENARIO_ID==5083){
                if( (fp_s[87]=fopen("./data/landuse_ann_1860_hist_1860_1860_1661-2299.flt","rb"))==NULL ){
                    printf("No landuse_ann_1860_hist_1860_1860_1661-2299.flt\n");
                    exit(1);
                }else{
                    Flag_FOPEN[87] ++;
                }
            }

            if(SCENARIO_ID==5090 || SCENARIO_ID==5091){
                if( (fp_s[87]=fopen("./data/landuse_ann_1860_1860_1860_1860_1661-2299.flt","rb"))==NULL ){
                    printf("No landuse_ann_1860_1860_1860_1860_1661-2299.flt\n");
                    exit(1);
                }else{
                    Flag_FOPEN[87] ++;
                }
            }else if(SCENARIO_ID==5041 || SCENARIO_ID==5047){
                if( (fp_s[87]=fopen("./data/landuse_ann_1860_hist_2005_2005_1661-2299.flt","rb"))==NULL ){
                    printf("No landuse_ann_1860_hist_2005_2005_1661-2299.flt\n");
                    exit(1);
                }else{
                    Flag_FOPEN[87] ++;
                }
            }else if(SCENARIO_ID==5043 || SCENARIO_ID==5045){
                if( (fp_s[87]=fopen("./data/landuse_ann_1860_hist_rcp26_2100rcp26_1661-2299_hadgem.flt","rb"))==NULL ){
                    printf("No landuse_ann_1860_hist_rcp26_2100rcp26_1661-2299_hadgem.flt\n");
                    exit(1);
                }else{
                    Flag_FOPEN[87] ++;
                }
            }else if(SCENARIO_ID==5042 || SCENARIO_ID==5092){
                if( (fp_s[87]=fopen("./data/landuse_ann_1860_hist_2005_2005_1661-2299.flt","rb"))==NULL ){
                    printf("No landuse_ann_1860_hist_2005_2005_1661-2299.flt\n");
                    exit(1);
                }else{
                    Flag_FOPEN[87] ++;
                }
            }else if(SCENARIO_ID==5044 || SCENARIO_ID==5046){  /* revised 2018/03/22 */
                if( (fp_s[87]=fopen("./data/landuse_ann_1860_hist_rcp60_2100_1661-2299_hadgem.flt","rb"))==NULL ){
                    printf("No landuse_ann_1860_hist_2005_2005_rcp60_2100_hadgem.flt\n");
                    exit(1);
                }else{
                    Flag_FOPEN[87] ++;
                }
            }else if(SCENARIO_ID==5093){
                if( (fp_s[87]=fopen("./data/landuse_ann_1860_hist_1860_1860_1661-2299.flt","rb"))==NULL ){
                    printf("No landuse_ann_1860_hist_1860_1860_1661-2299.flt\n");
                    exit(1);
                }else{
                    Flag_FOPEN[87] ++;
                }
            }
        
        }else{
            printf("BAD land-use setting for ISI-MIP2b\n");
            exit(1);
        }
    }else if(EX_BECCS == 2){
        /* BECCS data for GCP-MgNET-IIASA workshop, 2017/02/20 */
        if( (fp_s[87]=fopen("./data/beccs_201702.txt","rt"))==NULL ){
            printf("No beccs_201702.txt\n");
            exit(1); 
        }else{
            Flag_FOPEN[87] ++;
        }
    }else if(EX_BECCS == 3){
        /* BECCS data for Sustainability Science, 2017/10/27 */
        if( (fp_s[87]=fopen("./data/beccs_201710.txt","rt"))==NULL ){
            printf("No beccs_201710.txt\n");
            exit(1);
        }else{
            Flag_FOPEN[87] ++;
        }
    }else{
        /* BioFuel data: 2015/08/21 added by A.Ito  ***********************/
        /* data supplied from Kinoshita-san */
        if(BIOFUEL_RUN == 0){
            if( (fp_s[87]=fopen("./data/s10gcp_ssp_crop_2010-2100.txt","rt"))==NULL ){
                printf("No s10gcp_ssp_crop_2010-2100.txt\n");
                exit(1); 
            }else{
                Flag_FOPEN[87] ++;
            }
        }else if(BIOFUEL_RUN == 1){
            if( (fp_s[87]=fopen("./data/s10gcp_ssp_biofuel_2020-2100_current.txt","rt"))==NULL ){
                printf("No s10gcp_ssp_biofuel_2020-2100_current.txt\n");
                exit(1); 
            }else{
                Flag_FOPEN[87] ++;
            }
        }else if(BIOFUEL_RUN == 2){
            if( (fp_s[87]=fopen("./data/s10gcp_ssp_biofuel_2020-2100_low.txt","rt"))==NULL ){
                printf("No s10gcp_ssp_biofuel_2020-2100_low.txt\n");
                exit(1); 
            }else{
                Flag_FOPEN[87] ++;
            }
        }else if(BIOFUEL_RUN == 3){
            if( (fp_s[87]=fopen("./data/s10gcp_ssp_biofuel_2020-2100_middle.txt","rt"))==NULL ){
                printf("No s10gcp_ssp_biofuel_2020-2100_middle.txt\n");
                exit(1); 
            }else{
                Flag_FOPEN[87] ++;
            }
        }else{
            if( (fp_s[87]=fopen("./data/s10gcp_ssp_crop_2010-2100.txt","rt"))==NULL ){
                printf("No s10gcp_ssp_crop_2010-2100.txt\n");
                exit(1); 
            }else{
                Flag_FOPEN[87] ++;
            }
        }
    }
    
    /* updated: 2016/10/20 by A.Ito */
    if(SCENARIO_ID>=5010 && SCENARIO_ID<=5100){
        /* ISI-MIP2b: 2016/12/24 by A.Ito */
        if( (fp_s[88]=fopen("./data/nfert_ann_soc_1661-2299.flt","rb"))==NULL ){
            printf("No nfert_ann_soc_1661-2299.flt\n");
            exit(1); 
        }else{
            Flag_FOPEN[88] ++;
        }
    }else{
        /* NMIP: nitrogen input, 2015/11/19 by A.Ito ************/
        if( (fp_s[88]=fopen("./data/fin_nmip_v3.txt","rt"))==NULL ){
            /* printf("No fin_nmip_v2.txt\n"); */
            printf("No fin_nmip_v3.txt\n"); /* updated: 2017/10/17 by A.Ito */
            exit(1);
        }else{
            Flag_FOPEN[88] ++;
        }
    }
    
    /* N fertilizer scenario: 2016/11/22 by A.Ito **/
    if(EX_NFERT == 1){
        if( (fp_s[89]=fopen("./data/estni_hadgem_rcp2p6_ssp1_v1.txt","rt"))==NULL ){
            printf("No estni_hadgem_rcp2p6_ssp1_v1.txt\n");
            exit(1); 
        }else{
            Flag_FOPEN[89] ++;
        }
    }else if(EX_NFERT == 2){
        if( (fp_s[89]=fopen("./data/estni_hadgem_rcp2p6_ssp2_v1.txt","rt"))==NULL ){
            printf("No estni_hadgem_rcp2p6_ssp2_v1.txt\n");
            exit(1); 
        }else{
            Flag_FOPEN[89] ++;
        }
    }else if(EX_NFERT == 3){
        if( (fp_s[89]=fopen("./data/estni_hadgem_rcp2p6_ssp3_v1.txt","rt"))==NULL ){
            printf("No estni_hadgem_rcp2p6_ssp3_v1.txt\n");
            exit(1); 
        }else{
            Flag_FOPEN[89] ++;
        }
    }else if(EX_NFERT == 4){
        if( (fp_s[89]=fopen("./data/estni_hadgem_rcp4p5_ssp1_v1.txt","rt"))==NULL ){
            printf("No estni_hadgem_rcp4p5_ssp1_v1.txt\n");
            exit(1); 
        }else{
            Flag_FOPEN[89] ++;
        }
    }else if(EX_NFERT == 5){
        if( (fp_s[89]=fopen("./data/estni_hadgem_rcp4p5_ssp2_v1.txt","rt"))==NULL ){
            printf("No estni_hadgem_rcp4p5_ssp2_v1.txt\n");
            exit(1); 
        }else{
            Flag_FOPEN[89] ++;
        }
    }else if(EX_NFERT == 6){
        if( (fp_s[89]=fopen("./data/estni_hadgem_rcp4p5_ssp3_v1.txt","rt"))==NULL ){
            printf("No estni_hadgem_rcp4p5_ssp3_v1.txt\n");
            exit(1); 
        }else{
            Flag_FOPEN[89] ++;
        }
    }else if(EX_NFERT == 7){
        if( (fp_s[89]=fopen("./data/estni_miroc_rcp2p6_ssp1_v1.txt","rt"))==NULL ){
            printf("No estni_miroc_rcp2p6_ssp1_v1.txt\n");
            exit(1); 
        }else{
            Flag_FOPEN[89] ++;
        }
    }else if(EX_NFERT == 8){
        if( (fp_s[89]=fopen("./data/estni_miroc_rcp2p6_ssp2_v1.txt","rt"))==NULL ){
            printf("No estni_miroc_rcp2p6_ssp2_v1.txt\n");
            exit(1); 
        }else{
            Flag_FOPEN[89] ++;
        }
    }else if(EX_NFERT == 9){
        if( (fp_s[89]=fopen("./data/estni_miroc_rcp2p6_ssp3_v1.txt","rt"))==NULL ){
            printf("No estni_miroc_rcp2p6_ssp3_v1.txt\n");
            exit(1); 
        }else{
            Flag_FOPEN[89] ++;
        }
    }else if(EX_NFERT == 10){
        if( (fp_s[89]=fopen("./data/estni_miroc_rcp4p5_ssp1_v1.txt","rt"))==NULL ){
            printf("No estni_miroc_rcp4p5_ssp1_v1.txt\n");
            exit(1); 
        }else{
            Flag_FOPEN[89] ++;
        }
    }else if(EX_NFERT == 11){
        if( (fp_s[89]=fopen("./data/estni_miroc_rcp4p5_ssp2_v1.txt","rt"))==NULL ){
            printf("No estni_miroc_rcp4p5_ssp2_v1.txt\n");
            exit(1); 
        }else{
            Flag_FOPEN[89] ++;
        }
    }else if(EX_NFERT == 12){
        if( (fp_s[89]=fopen("./data/estni_miroc_rcp4p5_ssp3_v1.txt","rt"))==NULL ){
            printf("No estni_miroc_rcp4p5_ssp3_v1.txt\n");
            exit(1); 
        }else{
            Flag_FOPEN[89] ++;
        }
    }if(EX_NFERT == 101){
        /* Nishina ESSD data: nitrogen input, 2017/02/13 by A.Ito ********/
        if( (fp_s[89]=fopen("./data/ninput_201702.txt","rt"))==NULL ){
            printf("No ninput_201702.txt\n");
            exit(1); 
        }else{
            Flag_FOPEN[89] ++;
        }
    }else{
        if( (fp_s[89]=fopen("./data/estni_hadgem_rcp2p6_ssp1_v1.txt","rt"))==NULL ){
            printf("No estni_hadgem_rcp2p6_ssp1_v1.txt\n");
            exit(1); 
        }else{
            Flag_FOPEN[89] ++;
        }
    }
    
    if( (fp_s[90]=fopen("./data/bnd_impress.txt","rt"))==NULL ){
        printf("No bnd_impress.txt\n");
        exit(1); 
    }else{
        Flag_FOPEN[90] ++;
    }

    if( (fp_s[91]=fopen("./data/nfert_potter.txt","rt"))==NULL ){
        printf("No nfert_potter.txt\n");
        exit(1); 
    }else{
        Flag_FOPEN[91] ++;
    }

    if( (fp_s[92]=fopen("./data/wet_05_glwd_meris.txt","rt"))==NULL ){
        printf("No wet_05_glwd_meris.txt\n");
        exit(1);
    }else{
        Flag_FOPEN[92] ++;
    }
}
