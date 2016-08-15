/*	VISIT:	Vegetation Integrative SImulation Tool						*/
/*  Old name:	Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & RIGC/JAMSTEC						*/
/*  Carbon cycle, erosion, biomass burning, land-use change,			*/
/*  CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version 1.0.0	cerated in September 07, 2007						*/

#include<stdio.h>
#include<math.h>
#include"structure.h"
#include"prototype.h"

/* N budget ******************************************************/
void n_fertilizer_in(
	struct Grid *grid, 
	struct Loct *loct
){
    long nyear;
	double fert_input;
	extern double MDN[12];
	
	/* National average: from FAOSTAT, 2002-2003 */
	/* kg N ha-1 yr-1 */
	switch(grid->country){
		case 0:		fert_input = 57.65;	break;	/*	Ocean */
		case 4:		fert_input = 57.65;	break;	/*Afghanistan */
		case 8:		fert_input = 42.88;	break;	/*	Albania */
		case 12:	fert_input = 2.75;	break;	/*	Algeria */
		case 20:	fert_input = 57.65;	break;	/*	Andorra */
		case 24:	fert_input = 0.66;	break;	/*	Angola */
		/* case 660:	fert_input = ;	break; */	/*	Anguilla */
		case 28:	fert_input = 57.65;	break;	/*	Antigua-and-Barbuda */
		case 32:	fert_input = 22.08;	break;	/*	Argentina */
		case 51:	fert_input = 23.43;	break;	/*	Armenia */
		case 533:	fert_input = 57.65;	break;	/*	Aruba */
		case 36:	fert_input = 19.91;	break;	/*	Australia */
		case 40:	fert_input = 121.31;	break;	/*	Austria */
		case 31:	fert_input = 7.52;	break;	/*	Azerbaijan */
		case 44:	fert_input = 57.65;	break;	/*	Bahamas */
		case 48:	fert_input = 57.65;	break;	/*	Bahrain */
		case 50:	fert_input = 129.18;	break;	/*	Bangladesh */
		case 52:	fert_input = 37.06;	break;	/*	Barbados */
		case 112:	fert_input = 52.79;	break;	/*	Belarus */
		case 56:	fert_input = 57.65;	break;	/*	Belgium */
		case 84:	fert_input = 110.75;	break;	/*	Belize */
		case 204:	fert_input = 4.41;	break;	/*	Benin */
		case 60:	fert_input = 57.65;	break;	/*	Bermuda */
		case 64:	fert_input = 57.65;	break;	/*	Bhutan */
		case 68:	fert_input = 1.88;	break;	/*	Bolivia */
		case 72:	fert_input = 57.65;	break;	/*	Botswana */
		case 76:	fert_input = 30.91;	break;	/*	Brazil */
		case 96:	fert_input = 57.65;	break;	/*	Brunei */
		case 100:	fert_input = 83.18;	break;	/*	Bulgaria */
		case 854:	fert_input = 2.47;	break;	/*	Burkina */
		case 108:	fert_input = 0.29;	break;	/*	Burundi */
		case 116:	fert_input = 1.69;	break;	/*	Cambodia */
		case 120:	fert_input = 2.33;	break;	/*	Cameroon */
		case 124:	fert_input = 31.59;	break;	/*	Canada */
		/* case C04:	fert_input = ;	break; */	/*	Canary-Islands */
		case 132:	fert_input = 57.65;	break;	/*	Cape-Verde */
		case 136:	fert_input = 57.65;	break;	/*	Cayman-Islands */
		case 140:	fert_input = 57.65;	break;	/*	Central-African-Repu */
		case 148:	fert_input = 57.65;	break;	/*	Chad */
		case 152:	fert_input = 101.3;	break;	/*	Chile */
		case 156:	fert_input = 186.54;	break;	/*	China */
		case 170:	fert_input = 97.04;	break;	/*	Colombia */
		case 174:	fert_input = 57.65;	break;	/*	Comoros */
		case 178:	fert_input = 57.65;	break;	/*	Congo */
		case 184:	fert_input = 57.65;	break;	/*	Cook-Islands */
		case 188:	fert_input = 97.61;	break;	/*	Costa-Rica */
		case 192:	fert_input = 10.62;	break;	/*	Cuba */
		case 196:	fert_input = 72.56;	break;	/*	Cyprus */
		case 200:	fert_input = 92.34;	break;	/*	Czechoslovakia */
		case 208:	fert_input = 82.59;	break;	/*	Denmark */
		case 262:	fert_input = 57.65;	break;	/*	Djibouti */
		case 212:	fert_input = 57.65;	break;	/*	Dominica */
		case 214:	fert_input = 57.65;	break;	/*	Dominican-Republic */
		case 218:	fert_input = 49.54;	break;	/*	Ecuador */
		case 818:	fert_input = 390.69;	break;	/*	Egypt */
		case 222:	fert_input = 36.89;	break;	/*	El-Salvador */
		case 226:	fert_input = 57.65;	break;	/*	Equatorial-Guinea */
		case 233:	fert_input = 111.35;	break;	/*	Estonia */
		case 230:	fert_input = 6.04;	break;	/*	Ethiopia */
		case 238:	fert_input = 57.65;	break;	/*	Falkland-Islands */
		case 234:	fert_input = 57.65;	break;	/*	Faroe-Islands */
		case 242:	fert_input = 22.03;	break;	/*	Fiji */
		case 246:	fert_input = 77.46;	break;	/*	Finland */
		case 250:	fert_input = 116.93;	break;	/*	France */
		case 254:	fert_input = 57.65;	break;	/*	French-Guiana */
		case 258:	fert_input = 57.65;	break;	/*	French-Polynesia */
		case 266:	fert_input = 0.75;	break;	/*	Gabon */
		case 270:	fert_input = 66.89;	break;	/*	Gambia */
		/* case B02:	fert_input = ;	break; */	/*	Gaza-Strip */
		case 268:	fert_input = 21.42;	break;	/*	Georgia */
		case 276:	fert_input = 150.46;	break;	/*	Germany */
		case 288:	fert_input = 1.77;	break;	/*	Ghana */
		case 292:	fert_input = 57.65;	break;	/*	Gibraltar */
		case 300:	fert_input = 72.52;	break;	/*	Greece */
		case 304:	fert_input = 57.65;	break;	/*	Greenland */
		case 308:	fert_input = 57.65;	break;	/*	Grenada */
		case 312:	fert_input = 57.65;	break;	/*	Guadeloupe */
		case 320:	fert_input = 49.58;	break;	/*	Guatemala */
		/* case B01:	fert_input = ;	break; */	/*	Guernsey */
		case 324:	fert_input = 1.00;	break;	/*	Guinea */
		case 624:	fert_input = 57.65;	break;	/*	Guinea-Bissau */
		case 328:	fert_input = 22.77;	break;	/*	Guyana */
		case 332:	fert_input = 57.65;	break;	/*	Haiti */
		case 340:	fert_input = 18.30;	break;	/*	Honduras */
		case 344:	fert_input = 57.65;	break;	/*	Hong-Kong */
		case 348:	fert_input = 70.85;	break;	/*	Hungary */
		/* case B10:	fert_input = ;	break; */	/*	ISRAELI-OCCUPIED-TER */
		case 352:	fert_input = 1060.71;	break;	/*	Iceland */
		case 356:	fert_input = 62.40;	break;	/*	India */
		case 360:	fert_input = 66.64;	break;	/*	Indonesia */
		case 364:	fert_input = 46.29;	break;	/*	Iran */
		case 368:	fert_input = 57.65;	break;	/*	Iraq */
		case 536:	fert_input = 57.65;	break;	/*	Iraq-Saudi-Arabia-Ne */
		case 372:	fert_input = 301.13;	break;	/*	Ireland */
		case 376:	fert_input = 83.94;	break;	/*	Israel */
		case 380:	fert_input = 69.07;	break;	/*	Italy */
		case 384:	fert_input = 3.62;	break;	/*	Ivory-Coast */
		case 388:	fert_input = 14.96;	break;	/*	Jamaica */
		/* case C07:	fert_input = ;	break; */	/*	Jan-Mayen */
		case 392:	fert_input = 118.40;	break;	/*	Japan */
		/* case B03:	fert_input = ;	break; */	/*	Jersey */
		case 400:	fert_input = 264.40;	break;	/*	Jordan */
		case 398:	fert_input = 3.00;	break;	/*	Kazakhstan */
		case 404:	fert_input = 17.45;	break;	/*	Kenya */
		case 296:	fert_input = 57.65;	break;	/*	Kiribati */
		case 414:	fert_input = 587.78;	break;	/*	Kuwait */
		case 417:	fert_input = 8.48;	break;	/*	Kyrgyzstan */
		case 418:	fert_input = 57.65;	break;	/*	Laos */
		case 428:	fert_input = 48.84;	break;	/*	Latvia */
		case 422:	fert_input = 69.08;	break;	/*	Lebanon */
		case 426:	fert_input = 57.65;	break;	/*	Lesotho */
		case 430:	fert_input = 57.65;	break;	/*	Liberia */
		case 434:	fert_input = 25.23;	break;	/*	Libya */
		case 438:	fert_input = 57.65;	break;	/*	Liechtenstein */
		case 440:	fert_input = 16.41;	break;	/*	Lithuania */
		case 442:	fert_input = 417.94;	break;	/*	Luxembourg */
		case 446:	fert_input = 57.65;	break;	/*	Macau */
		case 450:	fert_input = 0.95;	break;	/*	Madagascar */
		case 454:	fert_input = 14.62;	break;	/*	Malawi */
		case 458:	fert_input = 61.73;	break;	/*	Malaysia */
		case 462:	fert_input = 0.84;	break;	/*	Maldives */
		case 466:	fert_input = 57.65;	break;	/*	Mali */
		case 470:	fert_input = 53.54;	break;	/*	Malta */
		case 584:	fert_input = 57.65;	break;	/*	Marshall-Islands */
		case 474:	fert_input = 57.65;	break;	/*	Martinique */
		case 478:	fert_input = 57.65;	break;	/*	Mauritania */
		case 480:	fert_input = 92.85;	break;	/*	Mauritius */
		case 175:	fert_input = 57.65;	break;	/*	Mayotte */
		case 484:	fert_input = 43.78;	break;	/*	Mexico */
		case 583:	fert_input = 57.65;	break;	/*	Micronesia */
		case 498:	fert_input = 7.89;	break;	/*	Moldova */
		case 492:	fert_input = 57.65;	break;	/*	Monaco */
		case 496:	fert_input = 3.41;	break;	/*	Mongolia */
		case 500:	fert_input = 57.65;	break;	/*	Montserrat */
		case 504:	fert_input = 24.43;	break;	/*	Morocco */
		case 508:	fert_input = 4.73;	break;	/*	Mozambique */
		case 104:	fert_input = 0.54;	break;	/*	Myanmar */
		case 516:	fert_input = 2.27;	break;	/*	Namibia */
		case 524:	fert_input = 11.11;	break;	/*	Nepal */
		case 528:	fert_input = 263.17;	break;	/*	Netherlands */
		case 530:	fert_input = 57.65;	break;	/*	Netherlands-Antilles */
		case 540:	fert_input = 57.65;	break;	/*	New-Caledonia */
		case 554:	fert_input = 97.12;	break;	/*	New-Zealand */
		case 558:	fert_input = 16.64;	break;	/*	Nicaragua */
		case 562:	fert_input = 0.37;	break;	/*	Niger */
		case 566:	fert_input = 3.49;	break;	/*	Nigeria */
		case 570:	fert_input = 57.65;	break;	/*	Niue */
		case 408:	fert_input = 57.65;	break;	/*	North-Korea */
		case 580:	fert_input = 57.65;	break;	/*	Northern-Mariana-Isl */
		case 578:	fert_input = 116.79;	break;	/*	Norway */
		case 512:	fert_input = 57.65;	break;	/*	Oman */
		case 586:	fert_input = 109.31;	break;	/*	Pakistan */
		case 585:	fert_input = 57.65;	break;	/*	Palau-Islands */
		case 591:	fert_input = 23.59;	break;	/*	Panama */
		case 598:	fert_input = 57.65;	break;	/*	Papua-New-Guinea */
		case 600:	fert_input = 12.25;	break;	/*	Paraguay */
		case 604:	fert_input = 48.11;	break;	/*	Peru */
		case 608:	fert_input = 49.45;	break;	/*	Philippines */
		case 612:	fert_input = 57.65;	break;	/*	Pitcairn */
		case 616:	fert_input = 58.19;	break;	/*	Poland */
		case 620:	fert_input = 77.73;	break;	/*	Portugal */
		case 630:	fert_input = 57.65;	break;	/*	Puerto-Rico */
		case 634:	fert_input = 57.65;	break;	/*	Qatar */
		case 638:	fert_input = 57.65;	break;	/*	Reunion */
		case 642:	fert_input = 24.93;	break;	/*	Romania */
		case 643:	fert_input = 4.67;	break;	/*	Russia */
		case 646:	fert_input = 57.65;	break;	/*	Rwanda */
		case 674:	fert_input = 57.65;	break;	/*	San-Marino */
		case 678:	fert_input = 57.65;	break;	/*	Sao-Tome */
		case 682:	fert_input = 45.90;	break;	/*	Saudi-Arabia */
		case 686:	fert_input = 6.70;	break;	/*	Senegal */
		case 690:	fert_input = 57.65;	break;	/*	Seychelles */
		case 694:	fert_input = 57.65;	break;	/*	Sierra-Leone */
		case 702:	fert_input = 57.65;	break;	/*	Singapore */
		case 90:	fert_input = 57.65;	break;	/*	Solomon-Islands */
		case 706:	fert_input = 57.65;	break;	/*	Somalia */
		case 710:	fert_input = 27.37;	break;	/*	South-Africa */
		case 410:	fert_input = 124.53;	break;	/*	South-Korea */
		case 724:	fert_input = 57.80;	break;	/*	Spain */
		case 144:	fert_input = 88.46;	break;	/*	Sri-Lanka */
		/* case B07:	fert_input = ;	break; */	/*	St.-Helena */
		case 662:	fert_input = 57.65;	break;	/*	St.-Lucia */
		/* case C11:	fert_input = ;	break; */	/*	St.-Martin */
		case 670:	fert_input = 57.65;	break;	/*	St.-Vincent */
		case 736:	fert_input = 3.28;	break;	/*	Sudan */
		case 740:	fert_input = 73.34;	break;	/*	Suriname */
		case 744:	fert_input = 57.65;	break;	/*	Svalbard */
		case 748:	fert_input = 57.65;	break;	/*	Swaziland */
		case 752:	fert_input = 71.08;	break;	/*	Sweden */
		case 756:	fert_input = 94.17;	break;	/*	Switzerland */
		case 760:	fert_input = 40.53;	break;	/*	Syria */
		case 762:	fert_input = 57.65;	break;	/*	Tadzhikistan */
		case 158:	fert_input = 57.65;	break;	/*	Taiwan */
		case 834:	fert_input = 5.77;	break;	/*	Tanzania */
		case 764:	fert_input = 61.74;	break;	/*	Thailand */
		case 768:	fert_input = 2.10;	break;	/*	Togo */
		case 776:	fert_input = 57.65;	break;	/*	Tonga */
		case 780:	fert_input = 338.23;	break;	/*	Trinidad */
		case 788:	fert_input = 9.15;	break;	/*	Tunisia */
		case 792:	fert_input = 49.93;	break;	/*	Turkey */
		case 795:	fert_input = 57.65;	break;	/*	Turkmenistan */
		case 796:	fert_input = 57.65;	break;	/*	Turks-And-Caicos-Isl */
		case 800:	fert_input = 0.55;	break;	/*	Uganda */
		case 804:	fert_input = 11.56;	break;	/*	Ukraine */
		case 784:	fert_input = 57.65;	break;	/*	United-Arab-Emirates */
		case 826:	fert_input = 197.89;	break;	/*	United-Kingdom */
		case 840:	fert_input = 63.63;	break;	/*	United-States */
		case 858:	fert_input = 27.77;	break;	/*	Uruguay */
		case 860:	fert_input = 57.65;	break;	/*	Uzbekistan */
		case 548:	fert_input = 57.65;	break;	/*	Vanuatu */
		case 862:	fert_input = 85.99;	break;	/*	Venezuela */
		case 704:	fert_input = 134.48;	break;	/*	Vietnam */
		case 92:	fert_input = 57.65;	break;	/*	Virgin-Islands-(Brit */
		case 850:	fert_input = 57.65;	break;	/*	Virgin-Islands-(USA) */
		case 732:	fert_input = 57.65;	break;	/*	Western-Sahara */
		case 882:	fert_input = 57.65;	break;	/*	Western-Samoa */
		case 886:	fert_input = 4.86;	break;	/*	Yemen */
		case 890:	fert_input = 60.15;	break;	/*	Yugoslavia */
		case 180:	fert_input = 57.65;	break;	/*	Zaire */
		case 894:	fert_input = 57.65;	break;	/*	Zambia */
		case 716:	fert_input = 20.36;	break;	/*	Zimbabwe */
		case 2222:	fert_input = 57.65;	break;	/* Antarctica	 */
		default:	fert_input = 57.65;				
	}
	
	loct->n_frtlz_in = fert_input * MDN[grid->m] / 365.0;
    
    /* biofuel experiments: 2015/09/03 revised by A.Ito */
    if(BIOFUEL_RUN == 1){
        /* current */
        /* loct->n_frtlz_in = 80.0 * MDN[grid->m] / 365.0; */ /* 2015/08/27 by A.Ito */
        ;
    }
    if(BIOFUEL_RUN == 2){
        /* low */
        /* loct->n_frtlz_in = 40.0 * MDN[grid->m] / 365.0; */ /* 2015/08/27 by A.Ito */
        loct->n_frtlz_in = 10.0 * MDN[grid->m] / 365.0;
    }
    if(BIOFUEL_RUN == 3){
        /* medium */
        /* loct->n_frtlz_in = 120.0 * MDN[grid->m] / 365.0; */ /* 2015/08/27 by A.Ito */
        loct->n_frtlz_in = 60.0 * MDN[grid->m] / 365.0;
    }
    
    /* NMIP run: 2015/11/19 by A.Ito *****/
    if(NMIP_RUN >= 1){
        nyear = grid->niny;

        if(NMIP_RUN == 9){
            nyear = FDY_NINY+1;
        }
    
        if(grid->niny>=1900 && grid->niny<=2012){
            loct->n_frtlz_in = grid->nmip_nfert[nyear - 1900] * MDN[grid->m] / 365.0;
        }else if(grid->niny<1900){
            loct->n_frtlz_in = grid->nmip_nfert[1900 - 1900] * MDN[grid->m] / 365.0;
        }else if(grid->niny>2012){
            loct->n_frtlz_in = grid->nmip_nfert[2012 - 1900] * MDN[grid->m] / 365.0;
        }
    }
    
    /* 2016/07/25 by A.Ito */
    /* loct->n_frtlz_in *= 0.01; */
}
