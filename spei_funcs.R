#' Càlcul espacial de l'Evapotranspiració Potencial (PET) or Reference Evapotranspiration (ETo)
#'
#' @param tmed SpatRaster de temperatura mitjana mensual (Celsius). Necessari per 'thornthwaite'.
#' @param tmin SpatRaster de temperatura mínima mensual (Celsius). Necessari per 'hargreaves'.
#' @param tmax SpatRaster de temperatura màxima mensual (Celsius). Necessari per 'hargreaves'.
#' @param method Caràcter. Mètode de càlcul: "thornthwaite" o "hargreaves".
#' @param start_year Numèric. Any d'inici de la sèrie temporal.
#' @param cores Numèric. Nombre de nuclis per a la paral·lelització (per defecte 1).
#' @return Un SpatRaster amb la PET mensual.
#' @export
#' @import terra
#' @import SPEI
calc_spatial_pet_eto <- function(tmed = NULL, tmin = NULL, tmax = NULL, 
                             method = c("thornthwaite", "hargreaves"), 
                             start_year, cores = 1) {
  
  method <- match.arg(method)
  
  # 1. Determinació del raster base i extracció de latituds en graus (WGS84)
  base_rast <- if (method == "thornthwaite") tmed else tmin
  if (is.null(base_rast)) stop("Falten les dades de temperatura per al mètode seleccionat.")
  
  coords <- terra::crds(base_rast, df = FALSE, na.rm = FALSE)
  pts <- terra::vect(coords, type = "points", crs = terra::crs(base_rast))
  pts_ll <- terra::project(pts, "EPSG:4326")
  
  lats <- base_rast[[1]]
  terra::values(lats) <- terra::geom(pts_ll)[, "y"]
  names(lats) <- "Latitud"
  
  # 2. Definició de funcions internes per a terra::app
  if (method == "thornthwaite") {
    input_stack <- c(tmed, lats)
    
    pet_fun <- function(x, yr) {
      n_lyr <- length(x) - 1
      t_series <- x[1:n_lyr]
      lat_val <- x[n_lyr + 1]
      
      if (any(is.na(t_series)) || is.na(lat_val)) return(rep(NA, n_lyr))
      t_ts <- ts(t_series, frequency = 12, start = c(yr, 1))
      
      tryCatch({
        as.numeric(SPEI::thornthwaite(t_ts, lat = lat_val, na.rm = TRUE))
      }, error = function(e) rep(NA, n_lyr))
    }
    
  } else if (method == "hargreaves") {
    if (is.null(tmin) || is.null(tmax)) stop("Hargreaves requereix 'tmin' i 'tmax'.")
    # Apilem: Tmin(1..n), Tmax(1..n), Latitud(1)
    input_stack <- c(tmin, tmax, lats)
    
    pet_fun <- function(x, yr) {
      n_lyr <- (length(x) - 1) / 2
      tmin_s <- x[1:n_lyr]
      tmax_s <- x[(n_lyr + 1):(2 * n_lyr)]
      lat_val <- x[length(x)]
      
      if (any(is.na(tmin_s)) || any(is.na(tmax_s)) || is.na(lat_val)) return(rep(NA, n_lyr))
      
      tmin_ts <- ts(tmin_s, frequency = 12, start = c(yr, 1))
      tmax_ts <- ts(tmax_s, frequency = 12, start = c(yr, 1))
      
      tryCatch({
        as.numeric(SPEI::hargreaves(Tmin = tmin_ts, Tmax = tmax_ts, lat = lat_val, na.rm = TRUE))
      }, error = function(e) rep(NA, n_lyr))
    }
  }
  
  # 3. Execució matricial/paral·lela
  pet_raster <- terra::app(input_stack, fun = pet_fun, yr = start_year, cores = cores)
  names(pet_raster) <- names(base_rast)
  
  return(pet_raster)
}


#' Càlcul espacial de l'SPEI basat en cubs de dades (stars)
#'
#' @param prec SpatRaster de precipitació mensual.
#' @param pet_eto SpatRaster d'Evapotranspiració Potencial mensual.
#' @param scale Numèric. Escala temporal de l'SPEI.
#' @param start_year Numèric. Any d'inici de la sèrie temporal.
#' @param start_month Numèric. Mes d'inici de la sèrie temporal. Per defecte 1.
#' @param kernel Llista. Tipus de kernel i shift.
#' @param distribution Caràcter. Distribució probabilística.
#' @param fit Caràcter. Mètode d'ajust.
#' @return Un SpatRaster amb els valors de l'SPEI.
#' @export
calc_spatial_spei <- function(prec, pet_eto, scale = 3, 
                              start_year, start_month = 1,
                              kernel = list(type = 'rectangular', shift = 0), 
                              distribution = 'log-Logistic', 
                              fit = 'ub-pwm') {
  
  d_balance <- prec - pet_eto
  d_stars <- stars::st_as_stars(d_balance)
  
  spei_wrapper <- function(x) {
    if (all(is.na(x))) return(rep(NA, length(x)))
    
    # 1. Creem la TS obligatòriament perquè SPEI agrupi bé els mesos (gener amb gener...)
    x_ts <- ts(x, frequency = 12, start = c(start_year, start_month))
    
    res <- tryCatch({
      obj <- SPEI::spei(x_ts, scale = scale, kernel = kernel, 
                        distribution = distribution, fit = fit, na.rm = TRUE)
      as.numeric(obj$fitted)
    }, error = function(e) rep(NA, length(x)))
    
    return(res)
  }
  
  spei_stars <- stars::st_apply(d_stars, MARGIN = c("x", "y"), FUN = spei_wrapper)
  
  spei_raster <- as(spei_stars, "SpatRaster")
  names(spei_raster) <- names(prec)
  
  return(spei_raster)
}


