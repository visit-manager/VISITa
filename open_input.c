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
    if(ISIMIP_RUN==0){
        if(DL_CRU == 102){
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
        }else if(DL_CRU==106){
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
        }else if(DL_CRU==109){
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
        }else if(DL_CRU==111){
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
        }else if(DL_CRU==112){
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
        }else if(DL_CRU==113){
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
        }else if(DL_CRU==114){
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
        }else{
            printf("No CRU data\n");
            exit(1);
        }
    }else /* ISI-MIP: 2012/06/27 by A.Ito */
    /* ICARUS SSPs: 2016/08/14 by A.Ito */
    if(ISIMIP_RUN == 1){
        switch(GCM_ID){
            case 2001: case 2005:
                case 2201: case 2203: case 2208: case 2210:
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
                case 2202: case 2204: case 2206: case 2209: case 2211:
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
                case 2205: case 2207:
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
                case 2301: case 2303: case 2308: case 2310:
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
                case 2302: case 2304: case 2306: case 2309: case 2311:
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
                case 2305: case 2307:
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
                case 2101: case 2103: case 2108: case 2110:
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
                case 2102: case 2104: case 2106: case 2109: case 2111:
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
                case 2105: case 2107:
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
                case 2401: case 2403: case 2408: case 2410:
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
                case 2402: case 2404: case 2406: case 2409: case 2411:
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
                case 2405: case 2407:
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
                case 2501: case 2503: case 2508: case 2510:
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
                case 2502: case 2504: case 2506: case 2509: case 2511:
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
                case 2505: case 2507:
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
        switch(GCM_ID){
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
        switch(GCM_ID){
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
    }
    
    /*********************************************************************/
	if( (fp_s[0]=fopen("./data/nat_reg_05.dat","rt"))==NULL ){  
		printf("No nat_reg_05.dat\n");  
		exit(1); 
	}
	if( (fp_s[1]=fopen("./data/olson_biome3_05.dat","rt"))==NULL ){
		printf("No olson_biome3_05.dat\n");  
		exit(1); 
	}
	if( (fp_s[2]=fopen("./data/climate_05.dat","rt"))==NULL ){
		printf("No climate_05.dat\n");  
		exit(1); 
	}
	if( (fp_s[3]=fopen("./data/geo_prop_05.dat","rt"))==NULL ){
		printf("No geo_prop_05.dat\n");  
		exit(1); 
	}
	if( (fp_s[4]=fopen("./data/prec_cru.dat","rt"))==NULL ){
		printf("No prec_cru.dat\n");  
		exit(1); 
	}
	if( (fp_s[5]=fopen("./data/rivnum_trip_05.dat","rt"))==NULL ){
		printf("No rivnum_trip_05.dat\n");  
		exit(1); 
	}
	if( (fp_s[6]=fopen("./data/Albedo_soil_EX.dat","rt"))==NULL ){
		printf("No Albedo_soil_EX.dat\n");  
		exit(1); 
	}
	if( (fp_s[7]=fopen("./data/pcnt_orgmatter.dat","rt"))==NULL ){  
		printf("No pcnt_orgmatter.dat\n");  
		exit(1); 
	}
	if( (fp_s[8]=fopen("./data/fls_slope.dat","rt"))==NULL ){  
		printf("No fls_slope.dat\n");  
		exit(1); 
	}
	if( (fp_s[9]=fopen("./data/fk_erodibility.dat","rt"))==NULL ){  
		printf("No fk_erodibility.dat\n");  
		exit(1); 
	}
	if( (fp_s[10]=fopen("./data/potveg_sage_05.dat","rt"))==NULL ){  
		printf("No potveg_sage_05.dat\n");  
		exit(1); 
	}
	if( (fp_s[11]=fopen("./data/luc_1700-1992.dat","rt"))==NULL ){  
		printf("No luc_1700-1992.dat\n");  
		exit(1); 
	}
	if( (fp_s[12]=fopen("./data/FAO_soil.textur.1nnegl.dat","rt"))==NULL ){  
		printf("No FAO_soil.textur.1nnegl.dat\n");  exit(1); 
	}
	if( (fp_s[13]=fopen("./data/bulkdens.dat","rt"))==NULL ){  
		printf("No bulkdens.dat\n");  
		exit(1); 
	}
	if( (fp_s[14]=fopen("./data/fieldcap.dat","rt"))==NULL ){  
		printf("No fieldcap.dat\n");  
		exit(1); 
	}
	if( (fp_s[15]=fopen("./data/soil_field_cap0-150_30min.dat","rt"))==NULL ){  
		printf("No soil_field_cap0-150_30min.dat\n");  
		exit(1); 
	}
	if( (fp_s[16]=fopen("./data/soil_field_cap0-30_30min.dat","rt"))==NULL ){  
		printf("No soil_field_cap0-30_30min.dat\n");  
		exit(1); 
	}
	if( (fp_s[17]=fopen("./data/soil_bulk_dens0-30_30min.dat","rt"))==NULL ){  
		printf("No soil_bulk_dens0-30_30min.dat\n");  
		exit(1); 
	}
	if( (fp_s[18]=fopen("./data/soil_clay_perc0-30_30min.dat","rt"))==NULL ){  
		printf("No soil_clay_perc0-30_30min.dat\n");  
		exit(1); 
	}
	if( (fp_s[19]=fopen("./data/soil_sand_perc0-30_30min.dat","rt"))==NULL ){  
		printf("No soil_sand_perc0-30_30min.dat\n");  
		exit(1); 
	}
	if( (fp_s[20]=fopen("./data/rastert_soilph1.dat","rt"))==NULL ){  
		printf("No rastert_soilph1.dat\n");  
		exit(1); 
	}

	if(ALT_FWETLAND == 1){
		/* Alternative data (NASA/GISS): 2011/03/30 by A.Ito */
		if( (fp_s[21]=fopen("./data/fwetland_giss.txt","rt"))==NULL ){  
			printf("No fwetland_giss.txt\n");  
			exit(1); 
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
		}
	}
    /* WSL-Perogon welnad map: 2014/2/4 by A.Ito */
    if( (fp_s[83]=fopen("./data/fwetland_peregon.txt","rt"))==NULL ){
        printf("No fwetland_peregon.txt\n");  
        exit(1); 
    }
	
    /* GCP-CH4 inundation data: revised 2014/08/17 by A.Ito */
    if(ALT_INUND == 3){
        if( (fp_s[84]=fopen("./data/fw_00_swe_v1-2_smoothed_1999-2013.txt","rt"))==NULL ){
            printf("No GCP-CH4 inundation data\n");
            exit(1);
        }
    }else if(ALT_INUND == 4){
        if( (fp_s[84]=fopen("./data/fw_00_swe_v1-2_smoothed_1999-2013_norice.txt","rt"))==NULL ){
            printf("No GCP-CH4 inundation data\n");
            exit(1);
        }
    }else if(ALT_INUND == 5){
        if( (fp_s[84]=fopen("./data/fw_00_swe_v1-2_smoothed_1999-2013_norice-nowater.txt","rt"))==NULL ){
            printf("No GCP-CH4 inundation data\n");
            exit(1);
        }
    }else if(ALT_INUND == 6){
        if( (fp_s[84]=fopen("./data/fw_swamps-glwd_2000-2012.txt","rt"))==NULL ){
            printf("No fw_swamps-glwd_2000-2012.txt data\n");
            exit(1);
        }
    }else{
        if( (fp_s[84]=fopen("./data/fw_swamp-biascor_1999-2013.txt","rt"))==NULL ){
            printf("No GCP-CH4 inundation data\n");
            exit(1);
        }
	}

	if(ALT_FWETLAND == 1){
		/* Alternative data (SAGE): 2011/03/30 by A.Ito */
		if( (fp_s[22]=fopen("./data/fpaddy_iis_sage.txt","rt"))==NULL ){  
			printf("No major_crps_rice_0.50x0.50.asc\n");  
			exit(1); 
		}
	}else{
		if( (fp_s[22]=fopen("./data/rice_sage_05.dat","rt"))==NULL ){  
			printf("No rice_sage_05.dat\n");  
			exit(1); 
		}
	}
	
	if( (fp_s[24]=fopen("./data/soil_n_1m.dat","rt"))==NULL ){  
		printf("No soil_n_1m.dat\n");  
		exit(1); 
	}
	
	/* Land use change: Hurtt et al. (2006) */
	if(LANDUSE == 6){
		if( (fp_s[26]=fopen("./data/EOS_SAGEHYDE_1D_crop.dat","rt"))==NULL ){  
			printf("No EOS_SAGEHYDE_1D_crop.dat\n");  
			exit(1); 
		}
		if( (fp_s[27]=fopen("./data/EOS_SAGEHYDE_1D_past.dat","rt"))==NULL ){  
			printf("No EOS_SAGEHYDE_1D_past.dat\n");  
			exit(1); 
		}
		if( (fp_s[28]=fopen("./data/EOS_SAGEHYDE_1D_prim.dat","rt"))==NULL ){  
			printf("No EOS_SAGEHYDE_1D_prim.dat\n");  
			exit(1); 
		}
		if( (fp_s[29]=fopen("./data/EOS_SAGEHYDE_1D_secd.dat","rt"))==NULL ){  
			printf("No EOS_SAGEHYDE_1D_secd.dat\n");  
			exit(1); 
		}
		if( (fp_s[30]=fopen("./data/EOS_SAGEHYDE_1D_ssma.dat","rt"))==NULL ){  
			printf("No EOS_SAGEHYDE_1D_ssma.dat\n");  
			exit(1); 
		}
		if( (fp_s[31]=fopen("./data/EOS_SAGEHYDE_1D_ssmb.dat","rt"))==NULL ){  
			printf("No EOS_SAGEHYDE_1D_ssmb.dat\n");  
			exit(1); 
		}
		if( (fp_s[32]=fopen("./data/EOS_SAGEHYDE_1D_t_cp.dat","rt"))==NULL ){  
			printf("No EOS_SAGEHYDE_1D_t_cp.dat\n");  
			exit(1); 
		}
		if( (fp_s[33]=fopen("./data/EOS_SAGEHYDE_1D_t_cs.dat","rt"))==NULL ){  
			printf("No EOS_SAGEHYDE_1D_t_cs.dat\n");  
			exit(1); 
		}
		if( (fp_s[34]=fopen("./data/EOS_SAGEHYDE_1D_t_pc.dat","rt"))==NULL ){  
			printf("No EOS_SAGEHYDE_1D_t_pc.dat\n");  
			exit(1); 
		}
		if( (fp_s[35]=fopen("./data/EOS_SAGEHYDE_1D_t_ps.dat","rt"))==NULL ){  
			printf("No EOS_SAGEHYDE_1D_t_ps.dat\n");  
			exit(1); 
		}
		if( (fp_s[36]=fopen("./data/EOS_SAGEHYDE_1D_t_sc.dat","rt"))==NULL ){  
			printf("No EOS_SAGEHYDE_1D_t_sc.dat\n");  
			exit(1); 
		}
		if( (fp_s[37]=fopen("./data/EOS_SAGEHYDE_1D_t_sp.dat","rt"))==NULL ){  
			printf("No EOS_SAGEHYDE_1D_t_sp.dat\n");  
			exit(1); 
		}
		if( (fp_s[38]=fopen("./data/EOS_SAGEHYDE_1D_t_ss1.dat","rt"))==NULL ){  
			printf("No EOS_SAGEHYDE_1D_t_ss1.dat\n");  
			exit(1); 
		}
		if( (fp_s[39]=fopen("./data/EOS_SAGEHYDE_1D_t_ss2.dat","rt"))==NULL ){  
			printf("No EOS_SAGEHYDE_1D_t_ss2.dat\n");  
			exit(1); 
		}
		if( (fp_s[40]=fopen("./data/EOS_SAGEHYDE_1D_t_ss3.dat","rt"))==NULL ){  
			printf("No EOS_SAGEHYDE_1D_t_ss3.dat\n");  
			exit(1); 
		}
		if( (fp_s[41]=fopen("./data/EOS_SAGEHYDE_1D_t_vc.dat","rt"))==NULL ){  
			printf("No EOS_SAGEHYDE_1D_t_vc.dat\n");  
			exit(1); 
		}
		if( (fp_s[42]=fopen("./data/EOS_SAGEHYDE_1D_t_vp.dat","rt"))==NULL ){  
			printf("No EOS_SAGEHYDE_1D_t_vp.dat\n");  
			exit(1); 
		}
		if( (fp_s[43]=fopen("./data/EOS_SAGEHYDE_1D_t_vs1.dat","rt"))==NULL ){  
			printf("No EOS_SAGEHYDE_1D_t_vs1.dat\n");  
			exit(1); 
		}
		if( (fp_s[44]=fopen("./data/EOS_SAGEHYDE_1D_t_vs2.dat","rt"))==NULL ){  
			printf("No EOS_SAGEHYDE_1D_t_vs2.dat\n");  
			exit(1); 
		}
	}
	/* U.NH harmonized historical land-use data (2010/01/31: by A.Ito) */
	/* URL  http://luh.unh.edu/ */
	if(LANDUSE==8){
		if( (fp_s[26]=fopen("./data/luc_eos2_gcrop_1700-2005.dat","rt"))==NULL ){  
			printf("No luc_eos2_gcrop_1700-2005.dat\n");  
			exit(1); 
		}
		if( (fp_s[27]=fopen("./data/luc_eos2_gpast_1700-2005.dat","rt"))==NULL ){  
			printf("No luc_eos2_gpast_1700-2005.dat\n");  
			exit(1); 
		}
		if( (fp_s[28]=fopen("./data/luc_eos2_gothr_1700-2005.dat","rt"))==NULL ){  
			printf("No luc_eos2_gothr_1700-2005.dat\n");  
			exit(1); 
		}
		if( (fp_s[29]=fopen("./data/luc_eos2_gsecd_1700-2005.dat","rt"))==NULL ){  
			printf("No luc_eos2_gsecd_1700-2005.dat\n");  
			exit(1); 
		}
		if( (fp_s[30]=fopen("./data/luc_eos2_ssma_1700-2005.dat","rt"))==NULL ){  
			printf("No luc_eos2_ssma_1700-2005.dat\n");  
			exit(1); 
		}
		if( (fp_s[31]=fopen("./data/luc_eos2_ssmb_1700-2005.dat","rt"))==NULL ){  
			printf("No luc_eos2_ssmb_1700-2005.dat\n");  
			exit(1); 
		}
		if( (fp_s[32]=fopen("./data/luc_eos2_gflcp_1700-2005.dat","rt"))==NULL ){  
			printf("No luc_eos2_gflcp_1700-2005.dat\n");  
			exit(1); 
		}
		if( (fp_s[33]=fopen("./data/luc_eos2_gflcs_1700-2005.dat","rt"))==NULL ){  
			printf("No luc_eos2_gflcs_1700-2005.dat\n");  
			exit(1); 
		}
		if( (fp_s[34]=fopen("./data/luc_eos2_gflpc_1700-2005.dat","rt"))==NULL ){  
			printf("No luc_eos2_gflpc_1700-2005.dat\n");  
			exit(1); 
		}
		if( (fp_s[35]=fopen("./data/luc_eos2_gflps_1700-2005.dat","rt"))==NULL ){  
			printf("No luc_eos2_gflps_1700-2005.dat\n");  
			exit(1); 
		}
		if( (fp_s[36]=fopen("./data/luc_eos2_gflsc_1700-2005.dat","rt"))==NULL ){  
			printf("No luc_eos2_gflsc_1700-2005.dat\n");  
			exit(1); 
		}
		if( (fp_s[37]=fopen("./data/luc_eos2_gflsp_1700-2005.dat","rt"))==NULL ){  
			printf("No luc_eos2_gflsp_1700-2005.dat\n");  
			exit(1); 
		}
		if( (fp_s[38]=fopen("./data/luc_eos2_gfsh1_1700-2005.dat","rt"))==NULL ){  
			printf("No luc_eos2_gfsh1_1700-2005.dat\n");  
			exit(1); 
		}
		if( (fp_s[39]=fopen("./data/luc_eos2_gfsh2_1700-2005.dat","rt"))==NULL ){  
			printf("No luc_eos2_gfsh2_1700-2005.dat\n");  
			exit(1); 
		}
		if( (fp_s[40]=fopen("./data/luc_eos2_gfsh3_1700-2005.dat","rt"))==NULL ){  
			printf("No luc_eos2_gfsh3_1700-2005.dat\n");  
			exit(1); 
		}
		if( (fp_s[41]=fopen("./data/luc_eos2_gflvc_1700-2005.dat","rt"))==NULL ){  
			printf("No luc_eos2_gflvc_1700-2005.dat\n");  
			exit(1); 
		}
		if( (fp_s[42]=fopen("./data/luc_eos2_gflvp_1700-2005.dat","rt"))==NULL ){  
			printf("No luc_eos2_gflvp_1700-2005.dat\n");  
			exit(1); 
		}
		if( (fp_s[43]=fopen("./data/luc_eos2_gfvh1_1700-2005.dat","rt"))==NULL ){  
			printf("No luc_eos2_gfvh1_1700-2005.dat\n");  
			exit(1); 
		}
		if( (fp_s[44]=fopen("./data/luc_eos2_gfvh2_1700-2005.dat","rt"))==NULL ){  
			printf("No luc_eos2_gfvh2_1700-2005.dat\n");  
			exit(1); 
		}
	}
    
    /* U.NH harmonized historical land-use data (2013/12/19: by A.Ito) */
	/* URL  http://luh.unh.edu/ */
	if(LANDUSE==9 || LANDUSE==10 || LANDUSE==11 || LANDUSE==12 || LANDUSE==13
        || LANDUSE==14 || LANDUSE==15 || LANDUSE==16 || LANDUSE==17){
        if( (fp_s[26]=fopen("./data/luh_gcrop_1500_2005.txt","rt"))==NULL ){  
			printf("No luh_gcrop_1500_2005.txt\n");  
			exit(1); 
		}
		if( (fp_s[27]=fopen("./data/luh_gpast_1500_2005.txt","rt"))==NULL ){  
			printf("No luh_gpast_1500_2005.txt\n");  
			exit(1); 
		}
		if( (fp_s[28]=fopen("./data/luh_gothr_1500_2005.txt","rt"))==NULL ){  
			printf("No luh_gothr_1500_2005.txt\n");  
			exit(1); 
		}
		if( (fp_s[29]=fopen("./data/luh_gsecd_1500_2005.txt","rt"))==NULL ){  
			printf("No luh_gsecd_1500_2005.txt\n");  
			exit(1); 
		}
		if( (fp_s[30]=fopen("./data/luh_gssma_1500_2005.txt","rt"))==NULL ){  
			printf("No luh_gssma_1500_2005.txt\n");  
			exit(1); 
		}
		if( (fp_s[31]=fopen("./data/luh_gssmb_1500_2005.txt","rt"))==NULL ){  
			printf("No luh_gssmb_1500_2005.txt\n");  
			exit(1); 
		}
		if( (fp_s[32]=fopen("./data/luh_gflcp_1500_2005.txt","rt"))==NULL ){  
			printf("No luh_gflcp_1500_2005.txt\n");  
			exit(1); 
		}
		if( (fp_s[33]=fopen("./data/luh_gflcs_1500_2005.txt","rt"))==NULL ){  
			printf("No luh_gflcs_1500_2005.txt\n");  
			exit(1); 
		}
		if( (fp_s[34]=fopen("./data/luh_gflpc_1500_2005.txt","rt"))==NULL ){  
			printf("No luh_gflpc_1500_2005.txt\n");  
			exit(1); 
		}
		if( (fp_s[35]=fopen("./data/luh_gflps_1500_2005.txt","rt"))==NULL ){  
			printf("No luh_gflps_1500_2005.txt\n");  
			exit(1); 
		}
		if( (fp_s[36]=fopen("./data/luh_gflsc_1500_2005.txt","rt"))==NULL ){  
			printf("No luh_gflsc_1500_2005.txt\n");  
			exit(1); 
		}
		if( (fp_s[37]=fopen("./data/luh_gflsp_1500_2005.txt","rt"))==NULL ){  
			printf("No luh_gflsp_1500_2005.txt\n");  
			exit(1); 
		}
		if( (fp_s[38]=fopen("./data/luh_gfsh1_1500_2005.txt","rt"))==NULL ){  
			printf("No luh_gfsh1_1500_2005.txt\n");  
			exit(1); 
		}
		if( (fp_s[39]=fopen("./data/luh_gfsh2_1500_2005.txt","rt"))==NULL ){  
			printf("No luh_gfsh2_1500_2005.txt\n");  
			exit(1); 
		}
		if( (fp_s[40]=fopen("./data/luh_gfsh3_1500_2005.txt","rt"))==NULL ){  
			printf("No luh_gfsh3_1500_2005.txt\n");  
			exit(1); 
		}
		if( (fp_s[41]=fopen("./data/luh_gflvc_1500_2005.txt","rt"))==NULL ){  
			printf("No luh_gflvc_1500_2005.txt\n");  
			exit(1); 
		}
		if( (fp_s[42]=fopen("./data/luh_gflvp_1500_2005.txt","rt"))==NULL ){  
			printf("No luh_gflvp_1500_2005.txt\n");  
			exit(1); 
		}
		if( (fp_s[43]=fopen("./data/luh_gfvh1_1500_2005.txt","rt"))==NULL ){  
			printf("No luh_gfvh1_1500_2005.txt\n");  
			exit(1); 
		}
		if( (fp_s[44]=fopen("./data/luh_gfvh2_1500_2005.txt","rt"))==NULL ){  
			printf("No luh_gfvh2_1500_2005.txt\n");  
			exit(1); 
		}
    }
    
    if(LANDUSE==9 || LANDUSE==10 || LANDUSE==18){
        if( (fp_s[59]=fopen("./data/luh_gcrop_rcp45_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gcrop_rcp45_2005_2100.txt\n");  
            exit(1); 
        }
        if( (fp_s[60]=fopen("./data/luh_gpast_rcp45_2005_2100.txt","rt"))==NULL ){
            printf("No luh_gpast_rcp45_2005_2100.txt\n");  
            exit(1); 
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
    }
    
    if(LANDUSE == 11 || LANDUSE == 17){
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
    }
    
    if(LANDUSE == 12){
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
    }
    
    if(LANDUSE == 13){
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
	
	/***************************************************/
	/* 0: stable */
	/* 1: gradual rise (SRES A1 or A1B) */
	/* 2: gradual rise (SRES A1FI) */
	/* 3: gradual rise (SRES A1T) */
	/* 4: gradual rise (SRES A2) */
	/* 5: gradual rise (SRES B1) */
	/* 6: gradual rise (SRES B2) */
	/* 7: uniform rise */
	if(GCM_ID == 0 || (GCM_ID>6000 && GCM_ID<7000)){
		CO2S = 4;
		fp_s[23] = fopen("./data/image_a2_fcrop.dat","rt");
		fp_s[45] = fopen("./data/image_a2_fgrass.dat","rt");
	}else if(GCM_ID ==1 || GCM_ID == 41 || GCM_ID == 61){
		CO2S = 1;
		fp_s[23]=fopen("./data/image_a1b_fcrop.dat","rt");
		fp_s[45]=fopen("./data/image_a1b_fgrass.dat","rt");
	}else if(GCM_ID==2){
		CO2S = 2;
		fp_s[23]=fopen("./data/image_a1b_fcrop.dat","rt");
		fp_s[45]=fopen("./data/image_a1b_fgrass.dat","rt");
	}else if(GCM_ID==3){
		CO2S = 3;
		fp_s[23]=fopen("./data/image_a1b_fcrop.dat","rt");
		fp_s[45]=fopen("./data/image_a1b_fgrass.dat","rt");
	}else if(GCM_ID==4 || GCM_ID==11 || GCM_ID==13 || GCM_ID==14 || GCM_ID==15 || GCM_ID==21 || GCM_ID==22 || 
			GCM_ID==23 || GCM_ID==31 || GCM_ID==42 || GCM_ID==51 || GCM_ID==62 || GCM_ID==71 ){
		CO2S = 4;
		fp_s[23]=fopen("./data/image_a2_fcrop.dat","rt");
		fp_s[45]=fopen("./data/image_a2_fgrass.dat","rt");
	}else if(GCM_ID==5 || GCM_ID==43){
		CO2S = 5;
		fp_s[23]=fopen("./data/image_b1_fcrop.dat","rt");
		fp_s[45]=fopen("./data/image_b1_fgrass.dat","rt");
	}else if(GCM_ID==6 || GCM_ID==12 || GCM_ID==16 || GCM_ID==17 || GCM_ID==18 || GCM_ID==24 || GCM_ID==32 || 
			GCM_ID==44 || GCM_ID==52 || GCM_ID==63 || GCM_ID==106){
		CO2S = 6;
		fp_s[23]=fopen("./data/image_b1_fcrop.dat","rt");
		fp_s[45]=fopen("./data/image_b1_fgrass.dat","rt");
	}
	
	if(GCM_ID==1000 ||GCM_ID==1010 ||GCM_ID==1011 ||GCM_ID==1012 ||GCM_ID==1050 ||GCM_ID==1060 ||GCM_ID==1070 ||
		GCM_ID==1080 ||GCM_ID==1090 ||GCM_ID==1100 ||GCM_ID==1101 ||GCM_ID==1110 ||GCM_ID==1120 ||
		GCM_ID==1130 ||GCM_ID==1131 ||GCM_ID==1132 ||GCM_ID==1133 ||GCM_ID==1134 ||GCM_ID==1150 ||
		GCM_ID==1151 ||GCM_ID==1152 ||GCM_ID==1153 ||GCM_ID==1160 ||GCM_ID==1161 ||GCM_ID==1162 ||
		GCM_ID==1170 ||GCM_ID==1180 ||GCM_ID==1190 ||GCM_ID==1200 ||GCM_ID==1201 ||GCM_ID==1202 ||
		GCM_ID==1203 ||GCM_ID==1204 ||GCM_ID==1220 ||GCM_ID==1221 ||GCM_ID==1222 ||GCM_ID==1223 ||
		GCM_ID==1224 ||GCM_ID==1225 ||GCM_ID==1226 ||GCM_ID==1240 ||GCM_ID==1250 ||GCM_ID==1251 ||
		GCM_ID==1260 ||GCM_ID==1261 ||GCM_ID==1262 ||GCM_ID==1263 ||GCM_ID==1270){ 
		CO2S = 1;
		
		if((fp_s[23]=fopen("./data/image_a1b_fcrop.dat","rt"))==NULL){
			printf("NO DATA:image_a1b_fcrop.dat !!\n");
			exit(1);
		}
        if((fp_s[45]=fopen("./data/image_a1b_fgrass.dat","rt"))==NULL){
			printf("NO DATA:image_a1b_fgrass.dat !!\n");
			exit(1);
        }
	}else if(GCM_ID==1013 ||GCM_ID==1014 ||GCM_ID==1015 || GCM_ID==1051 ||GCM_ID==1061 ||GCM_ID==1071 ||
		GCM_ID==1081 ||GCM_ID==1091 ||GCM_ID==1111 ||GCM_ID==1121 ||GCM_ID==1135 ||GCM_ID==1136 ||GCM_ID==1137 ||
		GCM_ID==1138 ||GCM_ID==1139 ||GCM_ID==1154 ||GCM_ID==1155 ||GCM_ID==1156 ||GCM_ID==1171 ||GCM_ID==1181 ||
		GCM_ID==1191 ||GCM_ID==1205 ||GCM_ID==1206 ||GCM_ID==1207 ||GCM_ID==1208 ||GCM_ID==1209 ||GCM_ID==1227 ||
		GCM_ID==1228 ||GCM_ID==1229 ||GCM_ID==1230 ||GCM_ID==1231 ||GCM_ID==1241 ||GCM_ID==1252 ||GCM_ID==1264 ||
		GCM_ID==1265 ||GCM_ID==1266 ||GCM_ID==1267 ||GCM_ID==1271){ 
		CO2S = 4;

		if((fp_s[23]=fopen("./data/image_a2_fcrop.dat","rt"))==NULL){
			printf("NO DATA !!\n");
			exit(1);
		}
		fp_s[45]=fopen("./data/image_a2_fgrass.dat","rt");
	}else if(GCM_ID==1001 ||GCM_ID==1016 ||GCM_ID==1017 ||GCM_ID==1018 ||GCM_ID==1052 ||GCM_ID==1062 ||GCM_ID==1072 ||
		GCM_ID==1082 ||GCM_ID==1092 ||GCM_ID==1102 ||GCM_ID==1103 ||GCM_ID==1112 ||GCM_ID==1122 ||GCM_ID==1140 ||
		GCM_ID==1141 ||GCM_ID==1142 ||GCM_ID==1143 ||GCM_ID==1144 ||GCM_ID==1157 ||GCM_ID==1158 ||GCM_ID==1159 ||
		GCM_ID==1163 ||GCM_ID==1164 ||GCM_ID==1165 ||GCM_ID==1182 ||GCM_ID==1192 ||GCM_ID==1210 ||GCM_ID==1211 ||
		GCM_ID==1212 ||GCM_ID==1213 ||GCM_ID==1214 ||GCM_ID==1232 ||GCM_ID==1233 ||GCM_ID==1234 ||GCM_ID==1235 ||
		GCM_ID==1236 ||GCM_ID==1237 ||GCM_ID==1238 ||GCM_ID==1239 ||GCM_ID==1242 ||GCM_ID==1253 ||GCM_ID==1272){ 
		CO2S = 5;
		
		if((fp_s[23]=fopen("./data/image_b1_fcrop.dat","rt"))==NULL){
			printf("NO DATA !!\n");
			exit(1);
		}
		fp_s[45]=fopen("./data/image_b1_fgrass.dat","rt");
	}
    
    /* ISI-MIP: 2012/06/27 by A.Ito ***********************************/
    /* assumption: SRES data were used, instead of RCP data */
    if(GCM_ID==2001 || GCM_ID==2011 || GCM_ID==2021 || GCM_ID==2031 || GCM_ID==2041
        || GCM_ID==2005 || GCM_ID==2015 || GCM_ID==2025 || GCM_ID==2035 || GCM_ID==2045
        || GCM_ID==4023){
        CO2S = 1; /* RCP2.6 */
        /**/
        if((fp_s[23]=fopen("./data/image_b1_fcrop.dat","rt"))==NULL){
			printf("NO image_b1_fcrop.dat !!\n");
			exit(1);
		}
		fp_s[45]=fopen("./data/image_b1_fgrass.dat","rt");
    }else if(GCM_ID==2002 || GCM_ID==2012 || GCM_ID==2022 || GCM_ID==2032 || GCM_ID==2042
        || GCM_ID==2006 || GCM_ID==2016 || GCM_ID==2026 || GCM_ID==2036 || GCM_ID==2046
         || GCM_ID==4012 || GCM_ID==4022){
        CO2S = 4; /* RCP8.5 */
        /**/
        if((fp_s[23]=fopen("./data/image_a2_fcrop.dat","rt"))==NULL){
			printf("NO image_a2_fcrop.dat !!\n");
			exit(1);
		}
		fp_s[45]=fopen("./data/image_a2_fgrass.dat","rt");
    }else if(GCM_ID==2003 || GCM_ID==2013 || GCM_ID==2023 || GCM_ID==2033 || GCM_ID==2043
        || GCM_ID==2007 || GCM_ID==2017 || GCM_ID==2027 || GCM_ID==2037 || GCM_ID==2047
         || GCM_ID==4011|| GCM_ID==4021){
        CO2S = 2; /* RCP4.5 */
        /**/
        if((fp_s[23]=fopen("./data/image_a1b_fcrop.dat","rt"))==NULL){
			printf("NO image_a1b_fcrop.dat !!\n");
			exit(1);
		}
		fp_s[45]=fopen("./data/image_a1b_fgrass.dat","rt");
    }else if(GCM_ID==2004 || GCM_ID==2014 || GCM_ID==2024 || GCM_ID==2034 || GCM_ID==2044
        || GCM_ID==2008 || GCM_ID==2018 || GCM_ID==2028 || GCM_ID==2038 || GCM_ID==2048
        || GCM_ID==4024){
        CO2S = 3; /* RCP6.0 */
        /**/
        if((fp_s[23]=fopen("./data/image_a1b_fcrop.dat","rt"))==NULL){
			printf("NO image_a1b_fcrop.dat !!\n");
			exit(1);
		}
		fp_s[45]=fopen("./data/image_a1b_fgrass.dat","rt");
    }
    
    /* GEO-MIP */
    if(GCM_ID>=3000 && GCM_ID<=3999){
        CO2S = 2; /* RCP4.5 */
        /**/
        if((fp_s[23]=fopen("./data/image_a1b_fcrop.dat","rt"))==NULL){
			printf("NO image_a1b_fcrop.dat !!\n");
			exit(1);
		}
		fp_s[45]=fopen("./data/image_a1b_fgrass.dat","rt");
    }
    
    /* ISI-MIP 2 (historical) */
    if(GCM_ID>=5000 && GCM_ID<=5999){
        CO2S = 2; /* RCP4.5 */
        /**/
        if((fp_s[23]=fopen("./data/image_a1b_fcrop.dat","rt"))==NULL){
			printf("NO image_a1b_fcrop.dat !!\n");
			exit(1);
		}
		fp_s[45]=fopen("./data/image_a1b_fgrass.dat","rt");
    }
	
     /* ISI-MIP: 2012/06/27 by A.Ito ***********************************/
    /* assumption: SRES data were used, instead of RCP data */
    if(GCM_ID==2101 || GCM_ID==2103 || GCM_ID==2108 || GCM_ID==2110 ||
        GCM_ID==2201 || GCM_ID==2203 || GCM_ID==2208 || GCM_ID==2210 ||
        GCM_ID==2301 || GCM_ID==2303 || GCM_ID==2308 || GCM_ID==2310 ||
        GCM_ID==2401 || GCM_ID==2403 || GCM_ID==2408 || GCM_ID==2410 ||
        GCM_ID==2501 || GCM_ID==2503 || GCM_ID==2508 || GCM_ID==2510){
        CO2S = 1; /* RCP2.6 */
        /**/
        
        switch(GCM_ID){
            case 2101:
                if((fp_s[23]=fopen("./data/fcrop_gfdl_aim_ssp1_rcp26_v1.txt","rt"))==NULL){
                    printf("NO fcrop_gfdl_aim_ssp1_rcp26_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2103:
                if((fp_s[23]=fopen("./data/fcrop_gfdl_aim_ssp2_rcp26_v1.txt","rt"))==NULL){
                    printf("NO fcrop_gfdl_aim_ssp2_rcp26_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2108:
                if((fp_s[23]=fopen("./data/fcrop_gfdl_aim_ssp4_rcp26_v1.txt","rt"))==NULL){
                    printf("NO fcrop_gfdl_aim_ssp4_rcp26_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2110:
                if((fp_s[23]=fopen("./data/fcrop_gfdl_aim_ssp5_rcp26_v1.txt","rt"))==NULL){
                    printf("NO fcrop_gfdl_aim_ssp5_rcp26_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2201:
                if((fp_s[23]=fopen("./data/fcrop_hadg_aim_ssp1_rcp26_v1.txt","rt"))==NULL){
                    printf("NO fcrop_hadg_aim_ssp1_rcp26_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2203:
                if((fp_s[23]=fopen("./data/fcrop_hadg_aim_ssp2_rcp26_v1.txt","rt"))==NULL){
                    printf("NO fcrop_hadg_aim_ssp2_rcp26_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2208:
                if((fp_s[23]=fopen("./data/fcrop_hadg_aim_ssp4_rcp26_v1.txt","rt"))==NULL){
                    printf("NO fcrop_hadg_aim_ssp4_rcp26_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2210:
                if((fp_s[23]=fopen("./data/fcrop_hadg_aim_ssp5_rcp26_v1.txt","rt"))==NULL){
                    printf("NO fcrop_hadg_aim_ssp5_rcp26_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2301:
                if((fp_s[23]=fopen("./data/fcrop_ipsl_aim_ssp1_rcp26_v1.txt","rt"))==NULL){
                    printf("NO fcrop_ipsl_aim_ssp1_rcp26_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2303:
                if((fp_s[23]=fopen("./data/fcrop_ipsl_aim_ssp2_rcp26_v1.txt","rt"))==NULL){
                    printf("NO fcrop_ipsl_aim_ssp2_rcp26_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2308:
                if((fp_s[23]=fopen("./data/fcrop_ipsl_aim_ssp4_rcp26_v1.txt","rt"))==NULL){
                    printf("NO fcrop_ipsl_aim_ssp4_rcp26_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2310:
                if((fp_s[23]=fopen("./data/fcrop_ipsl_aim_ssp5_rcp26_v1.txt","rt"))==NULL){
                    printf("NO fcrop_ipsl_aim_ssp5_rcp26_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2401:
                if((fp_s[23]=fopen("./data/fcrop_miro_aim_ssp1_rcp26_v1.txt","rt"))==NULL){
                    printf("NO fcrop_miro_aim_ssp1_rcp26_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2403:
                if((fp_s[23]=fopen("./data/fcrop_miro_aim_ssp2_rcp26_v1.txt","rt"))==NULL){
                    printf("NO fcrop_miro_aim_ssp2_rcp26_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2408:
                if((fp_s[23]=fopen("./data/fcrop_miro_aim_ssp4_rcp26_v1.txt","rt"))==NULL){
                    printf("NO fcrop_miro_aim_ssp4_rcp26_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2410:
                if((fp_s[23]=fopen("./data/fcrop_miro_aim_ssp5_rcp26_v1.txt","rt"))==NULL){
                    printf("NO fcrop_miro_aim_ssp5_rcp26_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2501:
                if((fp_s[23]=fopen("./data/fcrop_nore_aim_ssp1_rcp26_v1.txt","rt"))==NULL){
                    printf("NO fcrop_nore_aim_ssp1_rcp26_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2503:
                if((fp_s[23]=fopen("./data/fcrop_nore_aim_ssp2_rcp26_v1.txt","rt"))==NULL){
                    printf("NO fcrop_nore_aim_ssp2_rcp26_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2508:
                if((fp_s[23]=fopen("./data/fcrop_nore_aim_ssp4_rcp26_v1.txt","rt"))==NULL){
                    printf("NO fcrop_nore_aim_ssp4_rcp26_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2510:
                if((fp_s[23]=fopen("./data/fcrop_nore_aim_ssp5_rcp26_v1.txt","rt"))==NULL){
                    printf("NO fcrop_nore_aim_ssp5_rcp26_v1.txt !!\n");
                    exit(1);
                }
                break;
            default:
                printf("NO fcrop data !!\n");
                exit(1);
        }
        
		fp_s[45]=fopen("./data/image_b1_fgrass.dat","rt");
    }else if(GCM_ID==2102 || GCM_ID==2104 || GCM_ID==2106 || GCM_ID==2109 || GCM_ID==2111 ||
        GCM_ID==2202 || GCM_ID==2204 || GCM_ID==2206 || GCM_ID==2209 || GCM_ID==2211 ||
        GCM_ID==2302 || GCM_ID==2304 || GCM_ID==2306 || GCM_ID==2309 || GCM_ID==2311 ||
        GCM_ID==2402 || GCM_ID==2404 || GCM_ID==2406 || GCM_ID==2409 || GCM_ID==2411 ||
        GCM_ID==2502 || GCM_ID==2504 || GCM_ID==2506 || GCM_ID==2509 || GCM_ID==2511 ){
        CO2S = 2; /* RCP4.5 */
        /**/
        
        switch(GCM_ID){
            case 2102:
                if((fp_s[23]=fopen("./data/fcrop_gfdl_aim_ssp1_rcp45_v1.txt","rt"))==NULL){
                    printf("NO fcrop_gfdl_aim_ssp1_rcp45_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2104:
                if((fp_s[23]=fopen("./data/fcrop_gfdl_aim_ssp2_rcp45_v1.txt","rt"))==NULL){
                    printf("NO fcrop_gfdl_aim_ssp2_rcp45_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2106:
                if((fp_s[23]=fopen("./data/fcrop_gfdl_aim_ssp3_rcp45_v1.txt","rt"))==NULL){
                    printf("NO fcrop_gfdl_aim_ssp3_rcp45_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2109:
                if((fp_s[23]=fopen("./data/fcrop_gfdl_aim_ssp4_rcp45_v1.txt","rt"))==NULL){
                    printf("NO fcrop_gfdl_aim_ssp4_rcp45_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2111:
                if((fp_s[23]=fopen("./data/fcrop_gfdl_aim_ssp5_rcp45_v1.txt","rt"))==NULL){
                    printf("NO fcrop_gfdl_aim_ssp5_rcp45_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2202:
                if((fp_s[23]=fopen("./data/fcrop_hadg_aim_ssp1_rcp45_v1.txt","rt"))==NULL){
                    printf("NO fcrop_hadg_aim_ssp1_rcp45_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2204:
                if((fp_s[23]=fopen("./data/fcrop_hadg_aim_ssp2_rcp45_v1.txt","rt"))==NULL){
                    printf("NO fcrop_hadg_aim_ssp2_rcp45_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2206:
                if((fp_s[23]=fopen("./data/fcrop_hadg_aim_ssp3_rcp45_v1.txt","rt"))==NULL){
                    printf("NO fcrop_hadg_aim_ssp3_rcp45_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2209:
                if((fp_s[23]=fopen("./data/fcrop_hadg_aim_ssp4_rcp45_v1.txt","rt"))==NULL){
                    printf("NO fcrop_hadg_aim_ssp4_rcp45_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2211:
                if((fp_s[23]=fopen("./data/fcrop_hadg_aim_ssp5_rcp45_v1.txt","rt"))==NULL){
                    printf("NO fcrop_hadg_aim_ssp5_rcp45_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2302:
                if((fp_s[23]=fopen("./data/fcrop_ipsl_aim_ssp1_rcp45_v1.txt","rt"))==NULL){
                    printf("NO fcrop_ipsl_aim_ssp1_rcp45_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2304:
                if((fp_s[23]=fopen("./data/fcrop_ipsl_aim_ssp2_rcp45_v1.txt","rt"))==NULL){
                    printf("NO fcrop_ipsl_aim_ssp2_rcp45_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2306:
                if((fp_s[23]=fopen("./data/fcrop_ipsl_aim_ssp3_rcp45_v1.txt","rt"))==NULL){
                    printf("NO fcrop_ipsl_aim_ssp3_rcp45_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2309:
                if((fp_s[23]=fopen("./data/fcrop_ipsl_aim_ssp4_rcp45_v1.txt","rt"))==NULL){
                    printf("NO fcrop_ipsl_aim_ssp4_rcp45_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2311:
                if((fp_s[23]=fopen("./data/fcrop_ipsl_aim_ssp5_rcp45_v1.txt","rt"))==NULL){
                    printf("NO fcrop_ipsl_aim_ssp5_rcp45_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2402:
                if((fp_s[23]=fopen("./data/fcrop_miro_aim_ssp1_rcp45_v1.txt","rt"))==NULL){
                    printf("NO fcrop_miro_aim_ssp1_rcp45_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2404:
                if((fp_s[23]=fopen("./data/fcrop_miro_aim_ssp2_rcp45_v1.txt","rt"))==NULL){
                    printf("NO fcrop_miro_aim_ssp2_rcp45_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2406:
                if((fp_s[23]=fopen("./data/fcrop_miro_aim_ssp3_rcp45_v1.txt","rt"))==NULL){
                    printf("NO fcrop_miro_aim_ssp3_rcp45_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2409:
                if((fp_s[23]=fopen("./data/fcrop_miro_aim_ssp4_rcp45_v1.txt","rt"))==NULL){
                    printf("NO fcrop_miro_aim_ssp4_rcp45_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2411:
                if((fp_s[23]=fopen("./data/fcrop_miro_aim_ssp5_rcp45_v1.txt","rt"))==NULL){
                    printf("NO fcrop_miro_aim_ssp5_rcp45_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2502:
                if((fp_s[23]=fopen("./data/fcrop_nore_aim_ssp1_rcp45_v1.txt","rt"))==NULL){
                    printf("NO fcrop_nore_aim_ssp1_rcp45_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2504:
                if((fp_s[23]=fopen("./data/fcrop_nore_aim_ssp2_rcp45_v1.txt","rt"))==NULL){
                    printf("NO fcrop_nore_aim_ssp2_rcp45_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2506:
                if((fp_s[23]=fopen("./data/fcrop_nore_aim_ssp3_rcp45_v1.txt","rt"))==NULL){
                    printf("NO fcrop_nore_aim_ssp3_rcp45_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2509:
                if((fp_s[23]=fopen("./data/fcrop_nore_aim_ssp4_rcp45_v1.txt","rt"))==NULL){
                    printf("NO fcrop_nore_aim_ssp4_rcp45_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2511:
                if((fp_s[23]=fopen("./data/fcrop_nore_aim_ssp5_rcp45_v1.txt","rt"))==NULL){
                    printf("NO fcrop_nore_aim_ssp5_rcp45_v1.txt !!\n");
                    exit(1);
                }
                break;
            default:
                printf("NO fcrop data !!\n");
                exit(1);
        }
        
		fp_s[45]=fopen("./data/image_a1b_fgrass.dat","rt");
    }else if(GCM_ID==2105 || GCM_ID==2107 ||
        GCM_ID==2205 || GCM_ID==2207 ||
        GCM_ID==2305 || GCM_ID==2307 ||
        GCM_ID==2405 || GCM_ID==2407 ||
        GCM_ID==2505 || GCM_ID==2507){
        CO2S = 3; /* RCP6.0 */
        /**/
        switch(GCM_ID){
            case 2105:
                if((fp_s[23]=fopen("./data/fcrop_gfdl_aim_ssp2_rcp60_v1.txt","rt"))==NULL){
                    printf("NO fcrop_gfdl_aim_ssp2_rcp60_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2107:
                if((fp_s[23]=fopen("./data/fcrop_gfdl_aim_ssp3_rcp60_v1.txt","rt"))==NULL){
                    printf("NO fcrop_gfdl_aim_ssp3_rcp60_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2205:
                if((fp_s[23]=fopen("./data/fcrop_hadg_aim_ssp2_rcp60_v1.txt","rt"))==NULL){
                    printf("NO fcrop_hadg_aim_ssp2_rcp60_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2207:
                if((fp_s[23]=fopen("./data/fcrop_hadg_aim_ssp3_rcp60_v1.txt","rt"))==NULL){
                    printf("NO fcrop_hadg_aim_ssp3_rcp60_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2305:
                if((fp_s[23]=fopen("./data/fcrop_ipsl_aim_ssp2_rcp60_v1.txt","rt"))==NULL){
                    printf("NO fcrop_ipsl_aim_ssp2_rcp60_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2307:
                if((fp_s[23]=fopen("./data/fcrop_ipsl_aim_ssp3_rcp60_v1.txt","rt"))==NULL){
                    printf("NO fcrop_ipsl_aim_ssp3_rcp60_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2405:
                if((fp_s[23]=fopen("./data/fcrop_miro_aim_ssp2_rcp60_v1.txt","rt"))==NULL){
                    printf("NO fcrop_miro_aim_ssp2_rcp60_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2407:
                if((fp_s[23]=fopen("./data/fcrop_miro_aim_ssp3_rcp60_v1.txt","rt"))==NULL){
                    printf("NO fcrop_miro_aim_ssp3_rcp60_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2505:
                if((fp_s[23]=fopen("./data/fcrop_nore_aim_ssp2_rcp60_v1.txt","rt"))==NULL){
                    printf("NO fcrop_nore_aim_ssp2_rcp60_v1.txt !!\n");
                    exit(1);
                }
                break;
            case 2507:
                if((fp_s[23]=fopen("./data/fcrop_nore_aim_ssp3_rcp60_v1.txt","rt"))==NULL){
                    printf("NO fcrop_nore_aim_ssp3_rcp60_v1.txt !!\n");
                    exit(1);
                }
                break;
            default:
                printf("NO fcrop data !!\n");
                exit(1);
        }
        
		fp_s[45]=fopen("./data/image_a1b_fgrass.dat","rt");
    }
    
    /* GEO-MIP ***************************************/
    if(GCM_ID>=3000 && GCM_ID<=3999){
        CO2S = 2; /* RCP4.5 */
        /**/
        if((fp_s[23]=fopen("./data/image_a1b_fcrop.dat","rt"))==NULL){
			printf("NO image_a1b_fcrop.dat !!\n");
			exit(1);
		}
		fp_s[45]=fopen("./data/image_a1b_fgrass.dat","rt");
    }
    
    /* ISI-MIP 2 (historical) ******************************/
    if(GCM_ID>=5000 && GCM_ID<=5999){
        CO2S = 2; /* RCP4.5 */
        /**/
        if((fp_s[23]=fopen("./data/image_a1b_fcrop.dat","rt"))==NULL){
			printf("NO image_a1b_fcrop.dat !!\n");
			exit(1);
		}
		fp_s[45]=fopen("./data/image_a1b_fgrass.dat","rt");
    }

   /**************************************************************/
    /* SAGE crop data */
	if( (fp_s[46]=fopen("./data/sage_crops.dat","rt"))==NULL ){  
		printf("No sage_crops.dat\n");  
		exit(1); 
	}
    
    /* radiation diffusion parameters */
	if( (fp_s[47]=fopen("./data/reg_raddif_hd.dat","rt"))==NULL ){  
		printf("No diffuse rad model data\n");  
		exit(1); 
	}
	
    /* SSM/I inundation data */
	if( (fp_s[48]=fopen("./data/ssmi_season.dat","rt"))==NULL ){  
		printf("No SSMI inundation data\n");  
		exit(1); 
	}
    
    /* permafrost data */
	if( (fp_s[49]=fopen("./data/permafrost_nsidc.dat","rt"))==NULL ){  
		printf("No permafrost data\n");  
		exit(1); 
	}
	
	/* crop and pasture fractions: 1700-2007 */
	/* Ramankutty & Kimball: added 2010/07/07 (A.Ito) */
	if( (fp_s[50]=fopen("./data/glcrop_1700-2007_0.5.dat","rt"))==NULL ){  
		printf("No glcrop_1700-2007_0.5.dat\n");  
		exit(1); 
	}
	
	if( (fp_s[51]=fopen("./data/glpast_1700-2007_0.5.dat","rt"))==NULL ){  
		printf("No glpast_1700-2007_0.5.dat\n");  
		exit(1); 
	}
	
    /* N deposition **************************************************/
	/* by Galloway & Dentener */
    if( (fp_s[25]=fopen("./data/N-deposition.dat","rt"))==NULL ){
		printf("No N-deposition.dat\n");  
		exit(1); 
	}
	/* by CHASER */
	if( (fp_s[52]=fopen("./data/chaser_ndepo_2001mon.txt","rt"))==NULL ){  
		printf("No CHASER N deposition data: chaser_ndepo_2001mon.txt\n");
		exit(1); 
	}
	if( (fp_s[85]=fopen("./data/CHASER_V4.0+HTAP2_preind-1851-1855.txt","rt"))==NULL ){
		printf("No CHASER N deposition data: CHASER_V4.0+HTAP2_preind-1851-1855.txt\n");
		exit(1); 
	}
	if( (fp_s[86]=fopen("./data/CHASER_V4.0+CCMI_REFC1SD_2006-2010.txt","rt"))==NULL ){
		printf("No CHASER N deposition data: CHASER_V4.0+CCMI_REFC1SD_2006-2010.txt\n");
		exit(1); 
	}
	
	/* wood harvest based on RCP-harmonized data: LUHa.v1 *************/
	/* added by A.Ito (2010/10/15) */
    if(LANDUSE==9 || LANDUSE==10 || LANDUSE==11 || LANDUSE==12 || LANDUSE==13
         || LANDUSE==17){
        if( (fp_s[53]=fopen("./data/luh_gvbh1_1500_2005.txt","rt"))==NULL ){  
            printf("No luh_gvbh1_1500_2005.txt\n");  
            exit(1); 
        }
        if( (fp_s[54]=fopen("./data/luh_gvbh2_1500_2005.txt","rt"))==NULL ){  
            printf("No luh_gvbh2_1500_2005.txt\n");  
            exit(1); 
        }
        if( (fp_s[55]=fopen("./data/luh_gsbh1_1500_2005.txt","rt"))==NULL ){  
            printf("No luh_gsbh1_1500_2005.txt\n");  
            exit(1); 
        }
        if( (fp_s[56]=fopen("./data/luh_gsbh2_1500_2005.txt","rt"))==NULL ){  
            printf("No luh_gsbh2_1500_2005.txt\n");  
            exit(1); 
        }
        if( (fp_s[57]=fopen("./data/luh_gsbh3_1500_2005.txt","rt"))==NULL ){  
            printf("No luh_gsbh3_1500_2005.txt\n");  
            exit(1); 
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
        }
    }
        
    if(LANDUSE==9 || LANDUSE==10){
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
    }else if(LANDUSE==11 || LANDUSE==17){
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
    }else if(LANDUSE==12){
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
    }else if(LANDUSE==13){
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
    }else if(LANDUSE==14){
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
    }else if(LANDUSE==15){
        /* ICARUS SSP data by T.Kinoshita: 2014/09/04 => 09/18 ******************/
        if( (fp_s[78]=fopen("./data/luc_wh_ssp2_2001-2100.txt","rt"))==NULL ){
            printf("No luc_wh_ssp2_2001-2100.txt\n");
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
    }else if(LANDUSE==16){
        /* ICARUS SSP data by T.Kinoshita: 2014/09/04 => 09/18 ******************/
        if( (fp_s[78]=fopen("./data/luc_wh_ssp3_2001-2100.txt","rt"))==NULL ){
            printf("No luc_wh_ssp3_2001-2100.txt\n");
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
    }else{
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
    }
    
    /* GlobAlbedo *******************************/
    if( (fp_s[58]=fopen("./data/GlobAlbedo_av.flt","rb"))==NULL ){
        printf("No GlobAlbedo_av.flt\n");  
        exit(1); 
    }
    
    /* BioFuel data: 2015/08/21 added by A.Ito  ***********************/
    /* data supplied from Kinoshita-san */
    if(BIOFUEL_RUN==0){
        if( (fp_s[87]=fopen("./data/s10gcp_ssp_crop_2010-2100.txt","rt"))==NULL ){
            printf("No s10gcp_ssp_crop_2010-2100.txt\n");
            exit(1); 
        }
    }else if(BIOFUEL_RUN==1){
        if( (fp_s[87]=fopen("./data/s10gcp_ssp_biofuel_2020-2100_current.txt","rt"))==NULL ){
            printf("No s10gcp_ssp_biofuel_2020-2100_current.txt\n");
            exit(1); 
        }
    }else if(BIOFUEL_RUN==2){
        if( (fp_s[87]=fopen("./data/s10gcp_ssp_biofuel_2020-2100_low.txt","rt"))==NULL ){
            printf("No s10gcp_ssp_biofuel_2020-2100_low.txt\n");
            exit(1); 
        }
    }else if(BIOFUEL_RUN==3){
        if( (fp_s[87]=fopen("./data/s10gcp_ssp_biofuel_2020-2100_middle.txt","rt"))==NULL ){
            printf("No s10gcp_ssp_biofuel_2020-2100_middle.txt\n");
            exit(1); 
        }
    }else{
        if( (fp_s[87]=fopen("./data/s10gcp_ssp_crop_2010-2100.txt","rt"))==NULL ){
            printf("No s10gcp_ssp_crop_2010-2100.txt\n");
            exit(1); 
        }
    }
    
    /* NMIP: nitrogen input, 2015/11/19 by A.Ito ************/
    if( (fp_s[88]=fopen("./data/fin_nmip.txt","rt"))==NULL ){
        printf("No fin_nmip.txt\n");
        exit(1); 
    }
}
