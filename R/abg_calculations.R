# abg_calculations.R

#' Calculate A-a O2 Gradient
#'
#' This function calculates the Alveolar-arterial (A-a) gradient.
#'
#' @param paco2 A numeric, the partial pressure of carbon dioxide in arterial
#'   blood
#' @param pao2 A numeric, the partial pressure of oxygen in arterial blood
#' @param fio2 A numeric, the Fraction of Inspired Oxygen; defaults to room air,
#'   21 percent if NA
#' @param temp A numeric, the patient temperature in degrees Celsius; defaults
#'   to normal temperature, 37 degrees if NA
#' @param elev A numeric, the elevation above sea level in meters; defaults to
#'   sea level, 0 meters if NA
#'
#' @return A numeric
#' @export
#'
#' @references 1. Kanber GJ, King FW, Eshchar YR, et. al. The alveolar-arterial
#'   oxygen gradient in young and elderly men during air and oxygen breathing.
#'   Am Rev Respir Dis. 1968 Mar;97(3):376-81.
#'
#'   2. Mellemgaard K. The alveolar-arterial oxygen difference: its size and
#'   components in normal man. Acta Physiol Scand. 1966 May;67(1):10-20.
#'
aa_gradient_test <- function(paco2, pao2, fio2, temp = 37, elev = 0) {
  # Aa DO2 = (FIO2 * (Patm - PH2O) - (PaCO2 / 0.8)) - PaO2
  # Patm = 760 * exp(Elevation / -7000); Houston elevation = 43 feet (13.106 m)
  # PH2O = 47 * exp((Temp - 37) / 18.4)
  # Check fio2 variable if na, and default to 21% if missing
    fio2 <-  ifelse( is.na(fio2),.21,fio2)
  # Check temp variable if na, and default to 37C if missing
    temp <-  ifelse( is.na(temp),37,temp)
  # Check elev variable if na, and default to sea level if missing
    elev <-  ifelse( is.na(elev),0,elev)
  
  patm <- 760 * exp(elev / -7000)
  ph2o <- 47 * exp((temp - 37) / 18.4)
  
  fio2 <- purrr::map_dbl(fio2, ~ dplyr::if_else(is.na(.x),.21,dplyr::if_else( .x > 1, .x / 100, .x)))
  
  (fio2 * (patm - ph2o) - (paco2 / 0.8)) - pao2
}
  
#' This function calculates the Henderson-Hasselbach equation.
#'
#' @param paco2 A numeric, the partial pressure of carbon dioxide in arterial
#'   blood
#' @param hco3 A numeric, the calculated concentration of bicarb in the arterial 
#'   blood
#' @param ph A numeric, the acid/base balance in arterial blood, usually between
#'   7.35 and 7.45
#'
#' @return A numeric
#' @export
#'
#' @references 1. Samuel, T. R., Ilanchezian, Rajagopalan, B. Application of modified Henderson Equation 
#'   in abg interpretation Int. J. Pharm. Sci. Rev. Res., 37(2), March – April 2016; Article No. 30, 
#'   Pages: 169-177  
#'       
  
henderson_hasselbach <- function(paco2,hco3,ph) {
  # Henderson Hasselbach assesses whether the vvalues of the abg are consistent 
  # within themselves.  This is a validation tool for determining if the abg
  # is inaccurate.  Formula:
  #         
  #            H+ = (24*(PaCO2))/HCO3
  #  
  # 
  if (
    is.numeric(paco2) == FALSE || 
    is.numeric(hco3) == FALSE || 
    is.numeric(ph) == FALSE) {
        stop("paco2, hco3, and ph require numerical values")
  }
  # Need to define pH conversion to H+ in nmol/L
  h <- 10^-ph*1000000000
  range_h <- c((h-1),(h+1))
  
  if ( findInterval( 24*paco2/hco3, range_h) == 1 ) 
    {
    return( "ABG values are consistent")
    }
  else
  {
    return("ABG values are not consistent")
  }
  
}

#' This function calculates the Henderson-Hasselbach equation.
#'
#' @param ph A numeric, the acid/base balance in arterial blood, usually between
#'   7.35 and 7.45
#' @param paco2 A numeric, the partial pressure of carbon dioxide in arterial
#'   blood
#'
#' @return A two part text variable
#' @export
#'
#' @references  
#'  
#'  

acid_base_balance <- function(ph,paco2,hco3) {
  # Define Acidosis based on pH range
  acidosis <- range(c(1,7.34999999))
  # Define Alkalosis based on pH range
  alkalemia <- range(c(7.450000000001,10))
  # Define Alkalosis based on pH range
  balanced <- range(c(7.35,7.45))
  # Define Acidosis based on pH range
  co2_elevated <- range(c(40,1000))
  # Define Alkalosis based on pH range
  co2_decreased <- range(c(0,39.99999999))
  # Define Alkalosis based on pH range
  co2_normal <- range(c(35,45))
  # Define Acidosis based on pH range
  hco3_elevated <- range(c(24,1000))
  # Define Alkalosis based on pH range
  hco3_decreased <- range(c(0,24.99999999))
  # Define Alkalosis based on pH range
  hco3_normal <- range(c(22,26))
  
  if ( findInterval(ph,balanced) == 1)
  {
    # Define Acidosis when within normal pH range
    ph_norm_acidosis <- range(c(1,7.39999999))
    # Define Acidosis when within normal pH range
    ph_norm_alkalosis <- range( c(7.4,10))
    
    if ( findInterval(paco2,co2_normal) == 1 ||
         findInterval(hco3,hco3_normal) == 1 )
    {
      return( "ABG within Normal Ranges")
    }
    if ( findInterval(ph,ph_norm_acidosis) == 1 ||
         findInterval(paco2,co2_elevated) == 1 ||
         findInterval(hco3,hco3_elevated) == 1 )
    {
      return( "Respiratory Acidosis with metabolic compensation")
    }
    if ( findInterval(ph,ph_norm_alkalosis) == 1 ||
         findInterval(paco2,co2_decreased) == 1 ||
         findInterval(hco3,hco3_decreased) == 1 )
    {
      return( "Respiratory Alkalosis with metabolic compensation")
    }
  }
  if ( findInterval(ph,acidosis) == 1)
  {
    if (findInterval(paco2,co2_elevated) == 1 ||
        findInterval(hco3,hco3_normal) || 1)
    {
      return( "Respiratory Acidosis without metabolic compensation")
    }
    else
    {
      return( "Respiratory Acidosis with metabolic compensation")
    }
  }
  if ( findInterval(ph,alkalemia) == 1)
  {
    if (findInterval(paco2,co2_decreased) == 1 ||
        findInterval(hco3,hco3_normal) || 1)
    {
      return( "Respiratory Alkalosis without metabolic compensation")
    }
    else
    {
      return( "Respiratory Alkalosis with metabolic compensation")
    }
  }
}



#' This function calculates the anion gap equation.
#'
#' @param na A numeric, the serum/plasma sodium level in arterial or venous blood, usually between
#'   135-145, in mmol/L
#' @param cl A numeric, the serum/plasma chloride level in arterial or venous blood
#'   blood, in either mmol/L or mg/dL
#' @param bg A numeric, the serum/plasma blood glucose level in arterial or venous blood, usually between
#'   60-180 mg/dl or 4-12 mmol/L
#' @param osm A numeric, either serum/plasma osmolality  from arterial or venous blood
#' @param bun A numeric, serum/plasma BUN (in mg/dL) from arterial or venous blood
#' @param urea A numeric, serum/plasma BUN urea (in mmol/L) from arterial or venous blood
#' @param alb A numeric, serum/plasma albumin (in g/dL) from arterial or venous blood
#' @param phos A numeric, serum/plasma phos (in mg/dL) from arterial or venous blood
#'
#' @return A multi-part data-frame 
#' @export
#'
#' @references  1. Oh, M.S., Carroll, H.J. The anion gap. New England Journal of Medicine
#'    1977; 297(15): 814-817.   https://doi.org/10.1056/nejm197710132971507 PMID 895822
#'    
#'    2. Criner, G.J. Metabolic disturbances of acid-base and electrolytes. Critical Care Study Guide: 
#'    Text and Review. 2nd ed.  Philadelphia, PA: Springer; 2010:696.
#'    
#'    3. Kaufman, D.A. Interpretation of arterial bloog gases (ABGs).  American Thoracic Society.
#'    https://www.thoracic.org/professionals/clinical-resources/critical-care/clinical-education/abgs.php 
#'  
#'  

anion_gap <- function(na,cl,bicarb,bg,osm,bun = 0,urea = 0,alb = 4.5,phos = 3.5) {
  
  bg <- purrr::map_dbl(bg, ~ dplyr::if_else(.x > 20, .x / 18, .x))
  alb <- purrr::map_dbl(alb, ~ dplyr::if_else(is.na(.x),4.5,dplyr::if_else( .x > 10, .x / 10, .x)))
  phos <- purrr::map_dbl(phos, ~ dplyr::if_else( is.na(.x),3.5,dplyr::if_else( .x > 1.5, .x / 18, .x)))
  # Check elev variable if na, and default to sea level if missing
  # bun <-  ifelse( is.na(bun),0,bun)
  # Check elev variable if na, and default to sea level if missing
  # urea <-  ifelse( is.na(urea),0,urea)
  ur <- if_else(urea == 0, ( bun / 2.8), urea)
  
  ag <- ( na - ( cl + bicarb))
  alb_corrected_ag <- (ag + (2.5 * (4 - alb)))
  d_gap <- ag - 12
  alb_corr_d_gap <- (alb_corrected_ag - 12)
  d_ratio <- d_gap / (24 - bicarb)
  alb_corr_d_ratio <- alb_corr_d_gap / (24 - bicarb)

  # ag <- osm - (2 * na) - bg - ur
  osm_gap <- osm - (2*na - bg - ur)
  return( data.frame(anion_gap = ag,delta_gap = d_gap, delta_ratio = d_ratio, albumin_corr_ag = alb_corrected_ag, 
                     albumin_corr_delta_gap = alb_corr_d_gap, albumin_corr_delta_ratio = alb_corr_d_ratio,
                     osmolality_gap = osm_gap))
}

#' This function evaluates the anion gap equation.
#'
#' @param na A numeric, the serum/plasma sodium level in arterial or venous blood, usually between
#'   135-145, in mmol/L
#' @param cl A numeric, the serum/plasma chloride level in arterial or venous blood
#'   blood, in either mmol/L or mg/dL
#' @param bg A numeric, the serum/plasma blood glucose level in arterial or venous blood, usually between
#'   60-180 mg/dl or 4-12 mmol/L
#' @param osm A numeric, either serum/plasma osmolality  from arterial or venous blood
#' @param bun A numeric, serum/plasma BUN (in mg/dL) from arterial or venous blood
#' @param urea A numeric, serum/plasma BUN urea (in mmol/L) from arterial or venous blood
#'
#' @return A two part text variable
#' @export
#'
#' @references  
#'  
#'  


