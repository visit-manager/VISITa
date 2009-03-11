/*	VISIT: Vegetation Integrative SImulator for Trace gases				*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & EAIMG/ECRP/FRSGC					*/
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

extern long GCM, CO2S, GCM_R, GCM_C; /* */

/* open input files **************************************************/
void open_input(
	FILE *fp_s[IFILEN], 
	FILE *fp_c[4]
){
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
	if( (fp_s[21]=fopen("./data/wetland_glwd_05.dat","rt"))==NULL ){  
		printf("No wetland_glwd_05.dat\n");  
		exit(1); 
	}
	if( (fp_s[22]=fopen("./data/rice_sage_05.dat","rt"))==NULL ){  
		printf("No rice_sage_05.dat\n");  
		exit(1); 
	}
	if( (fp_s[24]=fopen("./data/soil_n_1m.dat","rt"))==NULL ){  
		printf("No soil_n_1m.dat\n");  
		exit(1); 
	}
	if( (fp_s[25]=fopen("./data/N-deposition.dat","rt"))==NULL ){  
		printf("No N-deposition.dat\n");  
		exit(1); 
	}
	
	/* Land use change: Hurtt et al. (2006) */
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

	/***************************************************/
	/* 0: stable */
	/* 1: gradual rise (SRES A1 or A1B) */
	/* 2: gradual rise (SRES A1FI) */
	/* 3: gradual rise (SRES A1T) */
	/* 4: gradual rise (SRES A2) */
	/* 5: gradual rise (SRES B1) */
	/* 6: gradual rise (SRES B2) */
	/* 7: uniform rise */
	if(GCM==0){
		CO2S = 4;
		fp_s[23]=fopen("./data/image_a2_fcrop.dat","rt");
		fp_s[45]=fopen("./data/image_a2_fgrass.dat","rt");
	}else if(GCM==1 || GCM==41 || GCM==61){
		CO2S = 1;
		fp_s[23]=fopen("./data/image_a1b_fcrop.dat","rt");
		fp_s[45]=fopen("./data/image_a1b_fgrass.dat","rt");
	}else if(GCM==2){
		CO2S = 2;
		fp_s[23]=fopen("./data/image_a1b_fcrop.dat","rt");
		fp_s[45]=fopen("./data/image_a1b_fgrass.dat","rt");
	}else if(GCM==3){
		CO2S = 3;
		fp_s[23]=fopen("./data/image_a1b_fcrop.dat","rt");
		fp_s[45]=fopen("./data/image_a1b_fgrass.dat","rt");
	}else if(GCM==4 || GCM==11 || GCM==13 || GCM==14 || GCM==15 || GCM==21 || GCM==22 || 
			GCM==23 || GCM==31 || GCM==42 || GCM==51 || GCM==62 || GCM==71 ){
		CO2S = 4;
		fp_s[23]=fopen("./data/image_a2_fcrop.dat","rt");
		fp_s[45]=fopen("./data/image_a2_fgrass.dat","rt");
	}else if(GCM==5 || GCM==43){
		CO2S = 5;
		fp_s[23]=fopen("./data/image_b1_fcrop.dat","rt");
		fp_s[45]=fopen("./data/image_b1_fgrass.dat","rt");
	}else if(GCM==6 || GCM==12 || GCM==16 || GCM==17 || GCM==18 || GCM==24 || GCM==32 || 
			GCM==44 || GCM==52 || GCM==63 || GCM==106){
		CO2S = 6;
		fp_s[23]=fopen("./data/image_b1_fcrop.dat","rt");
		fp_s[45]=fopen("./data/image_b1_fgrass.dat","rt");
	}
	
	if(GCM==1000 ||GCM==1010 ||GCM==1011 ||GCM==1012 ||GCM==1050 ||GCM==1060 ||GCM==1070 ||
		GCM==1080 ||GCM==1090 ||GCM==1100 ||GCM==1101 ||GCM==1110 ||GCM==1120 ||
		GCM==1130 ||GCM==1131 ||GCM==1132 ||GCM==1133 ||GCM==1134 ||GCM==1150 ||
		GCM==1151 ||GCM==1152 ||GCM==1153 ||GCM==1160 ||GCM==1161 ||GCM==1162 ||
		GCM==1170 ||GCM==1180 ||GCM==1190 ||GCM==1200 ||GCM==1201 ||GCM==1202 ||
		GCM==1203 ||GCM==1204 ||GCM==1220 ||GCM==1221 ||GCM==1222 ||GCM==1223 ||
		GCM==1224 ||GCM==1225 ||GCM==1226 ||GCM==1240 ||GCM==1250 ||GCM==1251 ||
		GCM==1260 ||GCM==1261 ||GCM==1262 ||GCM==1263){ 
		CO2S = 1;
		
		if((fp_s[23]=fopen("./data/image_a1b_fcrop.dat","rt"))==NULL){
			printf("NO DATA !!\n");
			exit(1);
		}
		fp_s[45]=fopen("./data/image_a1b_fgrass.dat","rt");
	}else if(GCM==1013 ||GCM==1014 ||GCM==1015 || GCM==1051 ||GCM==1061 ||GCM==1071 ||
		GCM==1081 ||GCM==1091 ||GCM==1111 ||GCM==1121 ||GCM==1135 ||GCM==1136 ||GCM==1137 ||
		GCM==1138 ||GCM==1139 ||GCM==1154 ||GCM==1155 ||GCM==1156 ||GCM==1171 ||GCM==1181 ||
		GCM==1191 ||GCM==1205 ||GCM==1206 ||GCM==1207 ||GCM==1208 ||GCM==1209 ||GCM==1227 ||
		GCM==1228 ||GCM==1229 ||GCM==1230 ||GCM==1231 ||GCM==1241 ||GCM==1252 ||GCM==1264 ||
		GCM==1265 ||GCM==1266 ||GCM==1267){ 
		CO2S = 4;

		if((fp_s[23]=fopen("./data/image_a2_fcrop.dat","rt"))==NULL){
			printf("NO DATA !!\n");
			exit(1);
		}
		fp_s[45]=fopen("./data/image_a2_fgrass.dat","rt");
	}else if(GCM==1001 ||GCM==1016 ||GCM==1017 ||GCM==1018 ||GCM==1052 ||GCM==1062 ||GCM==1072 ||
		GCM==1082 ||GCM==1092 ||GCM==1102 ||GCM==1103 ||GCM==1112 ||GCM==1122 ||GCM==1140 ||
		GCM==1141 ||GCM==1142 ||GCM==1143 ||GCM==1144 ||GCM==1157 ||GCM==1158 ||GCM==1159 ||
		GCM==1163 ||GCM==1164 ||GCM==1165 ||GCM==1182 ||GCM==1192 ||GCM==1210 ||GCM==1211 ||
		GCM==1212 ||GCM==1213 ||GCM==1214 ||GCM==1232 ||GCM==1233 ||GCM==1234 ||GCM==1235 ||
		GCM==1236 ||GCM==1237 ||GCM==1238 ||GCM==1239 ||GCM==1242 ||GCM==1253){ 
		CO2S = 5;
		
		if((fp_s[23]=fopen("./data/image_b1_fcrop.dat","rt"))==NULL){
			printf("NO DATA !!\n");
			exit(1);
		}
		fp_s[45]=fopen("./data/image_b1_fgrass.dat","rt");
	}
	
	if( (fp_s[46]=fopen("./data/sage_crops.dat","rt"))==NULL ){  
		printf("No sage_crops.dat\n");  
		exit(1); 
	}

	if( (fp_s[47]=fopen("./data/reg_raddif_hd.dat","rt"))==NULL ){  
		printf("No diffuse rad model data\n");  
		exit(1); 
	}
}
