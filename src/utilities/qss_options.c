/*
 * File:        qss_options.c
 * Description: Implementation file for Options structure.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>

#include "qss_options.h"
#include "qss_general_util.h"
#include "qss_spatial_derivatives3d.h"

#define DSZ sizeof(QSSLIB_REAL)

/*======== Helper Functions for Options structure manipulation ========*/
/*!
 * allocateOptions() allocates memory for an Options structure
 *
 * Arguments: none
 *
 * Return value: a pointer to the allocated memory
 *
 */
Options *allocateOptions(void);




/*================= Options structure manipulation ==================*/

Options *allocateOptions(void)
{
  Options *options = (Options *)calloc(1,sizeof(Options));

  return  options;
}


Options *createOptionsDefault(void)
{
  Options *options;
  
  options = allocateOptions();
  
  /* Main part of the structure */

  sprintf(options->outfile,"out_file");
  sprintf(options->path,"./");

  options->dx = 0.04;
  options->tmax = 150.0;
  options->tplot = 0.5;

  sprintf(options->accuracy,"medium");
  options->accuracy_id = MEDIUM;

  options->a = 0.15;
  options->b = 0.04;

  options->save_data = 1;
  options->do_reinit = 1;  
  options->do_mask = 1;

  options->eps_stop = 0.014;
  options->narrow_band = 0;
  options->overlap = 0;
  
  options->init_step = 0;
  options->print_details = 1;
  options->tmax_r = 10*options->dx;
  options->subcell_fix = 1;
  options->dc = 0.1;
  options->cmax = 1.0/ (options->dx/2.0);
  options->amax = 10;
  options->amin = 0.05;
  
  options->vol_frac_max = 0.97;
  options->vol_frac_min = 0;
  
  options->mask_reinit = 1;
  options->stop_touch = 1;

  options->narrow_band = 0;
  
  options->eps_coefficient = 1.5;

  options->extrapol_id = 0;
  options->extrapol_func=signedLinearExtrapolationBCqss;
  options->space_deriv_func=GET_RHS_3D_FIFTH;
  
  options->checkpoint = 0;
  
  options->reservoir_inlet = 1;

  options->theta = 0;
  options->C = 1;

  options->check_connectivity = 1;
  options->err_check_zone = 20*options->dx;
  
  options->phi_w_ind = 0;
  options->phi_nw_ind = 0;
  
  options->use_var_theta = 0;
  
  options->order_space_accur = 3;
  options->order_time_accur = 3;
  
  options->use_satn_stop = 1;
  options->conserve_imbibe = 0;
  options->center_inlet = 0;
  
  return options;
}


Options *createOptionsFromInputFile(char *filename)
{
  FILE    *fp;
  char    line[100], c, word[100], *exist, *new;
  Options *options;

  int     i, n, found, len, tmp1, cmp;
  int     dy_set=0, dz_set=0, b_set= 0, cmax_set=0;
  double  tmp;

  char *file_base;					 

  int N_accur_menu = sizeof(Accuracy_settings_menu)/
                     sizeof(Accuracy_settings_item);

  options = createOptionsDefault();

  if( (fp = fopen(filename,"r")) == NULL )
  {
     printf("\nCouldnt open file %s",filename);
     return options;
  }

  while( fgets(line,80,fp) != NULL )
  {
    /* skip empty spaces at the beginning of the line */
    n = 0;
    while(isspace(line[n])) n++;

    c = tolower(line[n]);

     /* Main part of the structure */

    if ( c == 'a' )
    { /* could be 'a' or 'accuracy' */
      if ( tolower(line[n+1]) == 'c' )
      {
        sscanf(line+n,"accuracy %s",word);
    
        /* try to identify input with one of the accuracy menu items */
        found = 0;
        for(i = 0; i < N_accur_menu; i++)
        {
          cmp = strcmp(word,Accuracy_settings_menu[i].name);
          if( cmp == 0 )
          {
            found = 1;
            options->accuracy_id = (QSSLIB_SPATIAL_DERIVATIVE_ACCURACY_TYPE) i;
            sprintf(options->accuracy,"%s",Accuracy_settings_menu[i].name); 
          }   
        }

        if( !found )
        {
          printf("\nAccuracy type %s not found, set to default.",word);
        }
      } else if ( tolower(line[n+2]) == 'a')
      {
        /* amax */
        sscanf(line+n,"%*s %lf ",&tmp);
        options->amax = (QSSLIB_REAL)tmp;
        
      } else if ( tolower(line[n+2]) == 'i')
      {
        /* amin */
        sscanf(line+n,"%*s %lf ",&tmp);
        options->amin = (QSSLIB_REAL)tmp;
        
      } else
      {
        sscanf(line+n,"%*s %lf ",&tmp);
        options->a = (QSSLIB_REAL)tmp;
      }
    }
    else if( c == 'b' )
    {  /* 'b' */
      sscanf(line+n,"%*s %lf ",&tmp);
      options->b = (QSSLIB_REAL)tmp; 
    }
    else if( c == 'c')
    {    
      if((tolower(line[n+1]) == 'm'))
      {
        /* cmax */
        sscanf(line+n,"%*s %lf ",&tmp);
        options->cmax = (QSSLIB_REAL)tmp;
	    cmax_set = 1;
      }
	  else if (tolower(line[n+5]) == 'p')
	  {  /* 'checkpoint' */
	      sscanf(line+n,"%*s %d ",&tmp1);   
	           
          if( (tmp1 == 0) || (tmp1 == 1))
             options->checkpoint = tmp1;
          else
	        printf("\nIncorrect checkpoint option %d, set to default.\n",tmp1);   
	   }
	   else if (tolower(line[n+6]) == 'c')
	  {  /* 'check_connectivity' */
	      sscanf(line+n,"%*s %d ",&tmp1);   
	      options->check_connectivity = tmp1;
             
	   }
	   else if (tolower(line[n+1]) == 'e')
	  {  /* 'center_inlet' */
	      sscanf(line+n,"%*s %d ",&tmp1);   
	      options->center_inlet = tmp1;
             
	   }
	}
    else if( c == 'd' )
    { /* could be 'dx' or 'do_reinit' */
      if ( tolower(line[n+1]) == 'x' )
      {
        sscanf(line+n,"%*s %lf ",&tmp);
        options->dx = (QSSLIB_REAL)tmp;
      }
      else if( line[n+3] == 'r')
      {
        sscanf(line+n,"%*s %d ",&tmp1);
        if ( (tmp1 == 0) || (tmp1 == 1))
          options->do_reinit = tmp1;
        else
        {
          printf("\nIncorrect do_reinit option %d, set to default.\n",tmp1);
        }    
      }
     else if( line[n+3] == 'm')
     { /* 'do_mask' */
	  sscanf(line+n,"%*s %d ",&tmp1);
	  if( (tmp1 == 0) || (tmp1 == 1))
             options->do_mask = tmp1;
	  else
	  {
	    printf("\nIncorrect do_mask option %d, set to default.\n",tmp1);
	  }    
     }
     else if( line[n+3] == 'r')
     { /* 'do_reinit' */
	  sscanf(line+n,"%*s %d ",&tmp1);
	  if( (tmp1 == 0) || (tmp1 == 1))
             options->do_reinit = tmp1;
	  else
	  {
	    printf("\nIncorrect do_reinit option %d, set to default.\n",tmp1);
	  }    
     }
     else if( tolower(line[n+1]) == 'c' )
     { /* 'dc' */
          sscanf(line+n,"%*s %lf ",&tmp);
          options->dc = (QSSLIB_REAL)tmp;
     } 
    }
	else if( c == 'e')
	{
	   if ( (tolower(line[n+1]) == 'p') && (tolower(line[n+4]) == 'c' ) )
	   { /* 'eps_coefficient' */
	        sscanf(line+n,"%*s %lf ",&tmp);
	        if( tmp > 0)
                options->eps_coefficient = (QSSLIB_REAL)tmp;
	        else
	        {
	            printf("\nIncorrect eps_coefficient option %g, set to default.\n",tmp);
	        }
	    }
        else if(  (tolower(line[n+1]) == 'r') && (tolower(line[n+10]) == 'z' ) )
	    { /* 'err_check_zone' */
	        sscanf(line+n,"%*s %lf ",&tmp);
            options->err_check_zone = (QSSLIB_REAL)tmp;
        }      
        else if(  (tolower(line[n+1]) == 'p') && (tolower(line[n+4]) == 's' ) )
	    {/* 'eps_stop' */
	        sscanf(line+n,"%*s %lf ",&tmp);
            options->eps_stop = (QSSLIB_REAL)tmp;
        }  
        else if(  (tolower(line[n+1]) == 'x') && (tolower(line[n+2]) == 't'))
        {   /* 'extrapol_id' */
            sscanf(line+n,"%*s %d ",&tmp1);
       
            switch(tmp1)
            {
                case 0: options->extrapol_id = 0;	         
		            options->extrapol_func=signedLinearExtrapolationBCqss;
		            break;
	            case 1: options->extrapol_id = 1;	         
		            options->extrapol_func=signedLinearExtrapolationBCqss;
		            break;
	            case 2: options->extrapol_id = 2;	         
		            options->extrapol_func=signedLinearExtrapolationBCqss;
		            break;
	            case 3: options->extrapol_id = 3;	         
		            options->extrapol_func=signedLinearExtrapolationBCqss;
		            break;
	            default:
	                printf("\nIncorrect extrapol_id option %d, set to default.\n",
		                                                          tmp1);
	            options->extrapol_id = 0;
                options->extrapol_func=signedLinearExtrapolationBCqss;
		        break;
            }
    
        }
    }
    else if ( c == 'i')
    {
        /* init_step */
        sscanf(line+n,"%*s %d ",&tmp1);
        options->init_step = tmp1;
    }
    else if( c == 'm' && line[n+1] == 'a')
    {   /* mask_reinit */
       sscanf(line+n,"%*s %d ",&tmp1);
       if( (tmp1 == 0) || (tmp1 == 1))
          options->mask_reinit = tmp1;
       else
       {
	      printf("\nIncorrect mask_reinit option %d, set to default.\n",tmp1);
       }    
    } 
    else if( c == 'n' )
    {  /* narrow_band */
       sscanf(line+n,"%*s %d ",&tmp1);
       if( (tmp1 == 0) || (tmp1 == 1))
          options->narrow_band = tmp1;
       else
       {
	      printf("\nIncorrect narrow_band option %d, set to default.\n",tmp1);
       }    
    }
    else if( c == 'o' && line[n+1] == 'u')
    { /* 'outfile' */
      sscanf(line+n,"outfile %s",word);
      sprintf(options->outfile,"%s",word);

      /* set path to the one given in outfile name */
      exist = strrchr(options->outfile,'/');
      if( exist )
      {
        len = strlen(options->outfile);
        for(i = len-1; i >= 0; i--)
        {
          if( (options->outfile)[i] == '/' ) break;
        }
        new  = malloc(i+1);
        strncpy(new,options->outfile,i+1);
        sprintf(options->path,"%s",new);
      }
      else sprintf(options->path,"./");
    }
    else if( c == 'o' && line[n+1] == 'v')
    { /* 'overlap' */
       sscanf(line+n,"%*s %lf ",&tmp);        
       options->overlap = tmp;
    }
    else if( c == 'r' && line[n+1] == 'e' )
    {  /* 'reservoir_inlet' */
       sscanf(line+n,"%*s %d ",&tmp1);
       if( (tmp1 == 0) || (tmp1 == 1))
           options->reservoir_inlet = tmp1;
       else
       {
	        printf("\nIncorrect reservoir_inlet option %d, set to default.\n",
	                                                           tmp1);
	        options->reservoir_inlet =  0;
       }    
    } 
    else if( c == 's' && (tolower(line[n+1]) == 'a') )
    { /* 'save_data' */  
      sscanf(line+n,"%*s %d ",&tmp1);
      if ( (tmp1 == 0) || (tmp1 == 1))
        options->save_data = tmp1;
      else
      {
        printf("\nIncorrect save_data option %d, set to default.\n",tmp1);
      }
    }
    else if( c == 's' && (tolower(line[n+5]) == 't') )
    {  /* 'stop_touch' */
       sscanf(line+n,"%*s %d ",&tmp1);
       if( (tmp1 == 0) || (tmp1 == 1))
          options->stop_touch = tmp1;
       else
       {
	      printf("\nIncorrect stop_touch option %d, set to default.\n",tmp1);
       }    
    }
    else if( c == 's' && (tolower(line[n+1]) == 'u') )
    {  /* 'subcell_fix' */
       sscanf(line+n,"%*s %d ",&tmp1);
       if( (tmp1 == 0) || (tmp1 == 1))
          options->subcell_fix = tmp1;
       else
       {
	      printf("\nIncorrect subcell_fix option %d, set to default.\n",tmp1);
       }    
    }
    else if( c == 's' && (tolower(line[n+1]) == 'p') )
    {  /* 'space_accur' */
       sscanf(line+n,"%*s %d ",&tmp1);
       options->order_space_accur = tmp1;   
    }
    else if( c == 't')
    { 
        if(tolower(line[n+1]) == 'm' )
        { /* tmax */
            sscanf(line+n,"%*s %lf ",&tmp);
            options->tmax = (QSSLIB_REAL)tmp;
        }
        else if( tolower(line[n+1]) == 'p' )
        {  /* 'tplot' */
            sscanf(line+n,"%*s %lf ",&tmp);
            options->tplot = (QSSLIB_REAL)tmp;
        }
        else if( tolower(line[n+1]) == 'i')
        {   /* time_accur */
           sscanf(line+n,"%*s %d ",&tmp1);
           options->order_time_accur = tmp1; 
        }
        else if( tolower(line[n+1]) == 'h')
        {	
	        /* contact angle, input in degrees */
	        sscanf(line+n,"%*s %lf ",&tmp);
	        options->theta = (QSSLIB_REAL)tmp * 3.14159/180;
	    }
	}
    else if( c == 'v') 
    {  
        if ((tolower(line[n+9]) == 'm' ) && (tolower(line[n+10]) == 'a' ))
        {/* 'vol_frac_max' */
            sscanf(line+n,"%*s %lf ",&tmp);
            options->vol_frac_max = (QSSLIB_REAL)tmp;
        } else if ((tolower(line[n+9]) == 'm' ) && (tolower(line[n+10]) == 'i' ))
        {/* 'vol_frac_min' */
            sscanf(line+n,"%*s %lf ",&tmp);
            options->vol_frac_min = (QSSLIB_REAL)tmp;
        }
    }
    else if( c == 'u')
    {
        if ((tolower(line[n+4]) == 'v')) {
            /* 'use_var_theta - for mixed wet conditions' */
            sscanf(line+n,"%*s %d ",&tmp1);
            options->use_var_theta = tmp1;
        } else if ((tolower(line[n+4]) == 's')) {
            /* use_satn_stop */
            sscanf(line+n,"%*s %d ",&tmp1);
            options->use_satn_stop = tmp1;
        }
    }

    /* User additions */
    /* Watch for overlaps and words starting with the same letter */
    
    if( c == 'p' )
    {  /* 'print_details' */
       sscanf(line+n,"%*s %d ",&tmp1);
       if( (tmp1 == 0) || (tmp1 == 1))
           options->print_details = tmp1;
       else
       {
	  printf("\nIncorrect print_details option %d, set to default.\n",
	                                                           tmp1);
	  options->print_details =  1;
       }    
    }
    
    /* end User additions */
  }

  fclose(fp);
  return options;
}


void  printOptions(Options *options,FILE *fp_out)
{

 fprintf(fp_out, "  -----------------------------------------------------------------------------\n");
 fprintf(fp_out, "  Geometry options: \n");
 fprintf(fp_out, "  -----------------------------------------------------------------------------\n");
 fprintf(fp_out, "  dx                   %8g  [ grid spacing ]\n",options->dx);
 fprintf(fp_out, "  mask_reinit          %8d  [ reinitialize mask or not ]\n",options->mask_reinit);
 fprintf(fp_out, "  reservoir_inlet      %8d  [ impose reservoir at inlet or not ]\n",options->reservoir_inlet); 
  

 fprintf(fp_out, "  -----------------------------------------------------------------------------\n");
 fprintf(fp_out, "  Running options: \n");
 fprintf(fp_out, "  -----------------------------------------------------------------------------\n");
 fprintf(fp_out, "  a                    %8g  [ Initial normal velocity (or cap pressure) ]\n", options->a);
 fprintf(fp_out, "  b                    %8g  [ Curvature coefficient   ]\n",options->b);
 fprintf(fp_out, "  dc                   %8g  [ curvature increment/decrement ]\n",options->dc); 
 fprintf(fp_out,"   space accuracy       %8d [ accuracy options: low, medium, high, very_high ]\n",options->order_space_accur);                                                               
 fprintf(fp_out,"   time accuracy        %8d  [ accuracy options: low, medium, high, very_high ]\n",options->order_time_accur);                                                               
 fprintf(fp_out, "  theta                %8g  [ contact angle, in radians ]\n",options->theta);
 fprintf(fp_out, "  tplot                %8g  [ intervals for reinit, error ]\n",options->tplot);
 fprintf(fp_out, "  init_step            %8d  [ Initial step (if restarting a stopped simulation) ]\n",options->init_step);
 fprintf(fp_out, "  overlap              %8g  [ overlap between pore space and mask allowed ]\n",options->overlap);
 fprintf(fp_out, "  eps_coefficient      %8g  [ numerical width of smearing for Heaviside ]\n",options->eps_coefficient); 
 fprintf(fp_out, "  check_connectivity   %8d  [ check for trapped components or not ]\n",options->check_connectivity); 
 fprintf(fp_out, "  checkpoint           %8d  [ output checkpoint data or not ]\n",options->checkpoint); 
 fprintf(fp_out, "  err_check_zone       %8g  [ distance from zero ls to check error ]\n",options->err_check_zone); 
 fprintf(fp_out, "  do_reinit            %8d  [ reinitilize periodically (1) or not (0)]\n",options->do_reinit);
 fprintf(fp_out, "  do_mask              %8d  [ impose mask (1) or not (0)]\n",options->do_mask);
 fprintf(fp_out, "  print_details        %8d  [ print details (1) or not (0)   ]\n",options->print_details);                                                                                        
 fprintf(fp_out, "  center_inlet         %8d  [ Inlet at center or not \n", options->center_inlet);


 fprintf(fp_out, "  -----------------------------------------------------------------------------\n");
 fprintf(fp_out, "  Stopping options: \n");
 fprintf(fp_out, "  -----------------------------------------------------------------------------\n");                                                                  
 fprintf(fp_out, "  tmax                 %8g  [ max running time allowed ]\n",options->tmax);
 fprintf(fp_out, "  eps_stop             %8g  [ criteria for stopping ]\n",options->eps_stop);
 fprintf(fp_out, "  tmax_reinit          %8g  [ max reinit time allowed ]\n",options->tmax_r);  
 fprintf(fp_out, "  cmax                 %8g  [ maximum curvature ]\n",options->cmax);
 fprintf(fp_out, "  amax                 %8g  [ maximum a, for drainage ]\n",options->amax);
 fprintf(fp_out, "  amin                 %8g  [ minimum a, for imbibition ]\n",options->amin); 
 fprintf(fp_out, "  vol_frac_max         %8g  [ maximum nw fraction allowed ]\n",options->vol_frac_max);
 fprintf(fp_out, "  stop_touch           %8d  [ stop the simulation when opp. bdry. is touched ]\n",options->stop_touch);                                          
 fprintf(fp_out, "  save_data            %8d  [ save data (1) or not (0)       ]\n",options->save_data);  
                                                                                                                       
 fprintf(fp_out, "  outfile              %8s  [ output file name ]\n",options->outfile); /* output file name */
 
}


void destroyOptions(Options *options)
{
  if(options) free(options);
}
