#' Cálculo de Tamaño Muestral para Estudios de Prevalencia
#' 
#' @param prevalencia Prevalencia verdadera esperada
#' @param precision Precisión absoluta deseada (requerida si se calcula tamaño muestral)
#' @param n Tamaño muestral (requerido si se calcula precisión)
#' @param nivel_confianza Nivel de confianza (por defecto 0.95)
#' @param N Tamaño de la población (usar NA para población infinita)
#' @param sensibilidad Sensibilidad de la prueba (por defecto 1)
#' @param especificidad Especificidad de la prueba (por defecto 1)
#' @param prop_estratos Proporciones de población en cada estrato (por defecto NULL para muestreo no estratificado)
#' @param efecto_diseno Efecto de diseño debido a estratificación (por defecto 1 para muestreo no estratificado)
#' @param calcular Qué calcular: "muestra" o "precision" (por defecto "muestra")
#' @return Una lista con los resultados calculados y parámetros de entrada
#' @export
CalculadoraPrevalencia <- function(prevalencia, precision = NULL, n = NULL, 
                                 nivel_confianza = 0.95, N = NA,
                                 sensibilidad = 1, especificidad = 1, 
                                 prop_estratos = NULL, efecto_diseno = 1,
                                 calcular = "muestra") {
  
  if (calcular != "muestra" && calcular != "precision") {
    stop("'calcular' debe ser 'muestra' o 'precision'")
  }
  
  if (calcular == "muestra" && is.null(precision)) {
    stop("Se debe proporcionar la precisión al calcular tamaño muestral")
  }
  
  if (calcular == "precision" && is.null(n)) {
    stop("Se debe proporcionar n al calcular precisión")
  }
  
  z <- qnorm((1 + nivel_confianza) / 2)
  
  # Prevalencia aparente
  Pa <- prevalencia * sensibilidad + (1 - prevalencia) * (1 - especificidad)
  
  if (calcular == "muestra") {
    if (is.na(N)) {
      # Fórmula para población infinita
      n <- (z^2 * Pa * (1 - Pa)) / (precision^2 * (sensibilidad + especificidad - 1)^2)
    } else {
      # Fórmula para población finita
      n <- (N * z^2 * Pa * (1 - Pa)) / ((precision^2 * (N - 1) * (sensibilidad + especificidad - 1)^2) + (z^2 * Pa * (1 - Pa)))
    }
    
    # Aplicar efecto de diseño
    n <- n * efecto_diseno
    
    # Calcular tamaños de estratos
    tamanos_estratos <- NULL
    if (!is.null(prop_estratos)) {
      if (abs(sum(prop_estratos) - 1) > 1e-6) {
        warning("Las proporciones de estratos no suman 1. Normalizando.")
        prop_estratos <- prop_estratos / sum(prop_estratos)
      }
      tamanos_estratos <- ceiling(n * prop_estratos)
      n <- sum(tamanos_estratos)
    }
    
    resultados <- list(
      tamano_muestral = ceiling(n),
      tamanos_estratos = tamanos_estratos,
      precision = precision
    )
  } else {
    # Ajustar tamaño muestral efectivo por efecto de diseño
    n_efectivo <- n / efecto_diseno
    
    if (is.na(N)) {
      # Fórmula para población infinita
      precision <- z * sqrt((Pa * (1 - Pa)) / n_efectivo) / (sensibilidad + especificidad - 1)
    } else {
      # Fórmula para población finita
      fpc <- sqrt((N - n_efectivo) / (N - 1))
      precision <- z * sqrt((Pa * (1 - Pa)) / n_efectivo) * fpc / (sensibilidad + especificidad - 1)
    }
    
    # Calcular tamaños de estratos
    tamanos_estratos <- NULL
    if (!is.null(prop_estratos)) {
      if (abs(sum(prop_estratos) - 1) > 1e-6) {
        warning("Las proporciones de estratos no suman 1. Normalizando.")
        prop_estratos <- prop_estratos / sum(prop_estratos)
      }
      tamanos_estratos <- ceiling(n * prop_estratos)
    }
    
    resultados <- list(
      precision = precision,
      tamano_muestral = n,
      tamanos_estratos = tamanos_estratos
    )
  }
  
  # Agregar elementos comunes
  resultados <- c(resultados, list(
    prevalencia_verdadera = prevalencia,
    prevalencia_aparente = Pa,
    nivel_confianza = nivel_confianza,
    tamano_poblacion = N,
    sensibilidad = sensibilidad,
    especificidad = especificidad,
    efecto_diseno = efecto_diseno,
    prop_estratos = prop_estratos
  ))
  
  return(resultados)
}

#' Cálculo Logístico para Estudios Transversales
#'
#' @param n_final Tamaño muestral final calculado
#' @param tasa_rechazo Tasa de rechazo esperada
#' @param tasa_elegibilidad Tasa de elegibilidad esperada
#' @param sujetos_por_dia Número de sujetos que pueden ser procesados por día
#' @param dias_laborables_mes Número de días laborables por mes
#' @return Una lista con los cálculos logísticos del estudio
#' @export
logistica_estudio_transversal <- function(n_final, 
                                        tasa_rechazo, 
                                        tasa_elegibilidad, 
                                        sujetos_por_dia, 
                                        dias_laborables_mes) {
  
  # Validación de entradas
  if (tasa_rechazo < 0 || tasa_rechazo >= 1) 
    stop("La tasa de rechazo debe estar entre 0 y 1")
  if (tasa_elegibilidad <= 0 || tasa_elegibilidad > 1) 
    stop("La tasa de elegibilidad debe estar entre 0 y 1")
  if (sujetos_por_dia <= 0) 
    stop("El número de sujetos por día debe ser positivo")
  if (dias_laborables_mes <= 0) 
    stop("El número de días laborables por mes debe ser positivo")
  
  n_evaluar <- n_final / (1 - tasa_rechazo)
  n_invitar <- n_evaluar / tasa_elegibilidad
  total_dias <- n_invitar / sujetos_por_dia
  duracion_meses <- total_dias / dias_laborables_mes
  
  resultados <- list(
    muestra_final = n_final,
    muestra_evaluar = ceiling(n_evaluar),
    muestra_invitar = ceiling(n_invitar),
    dias_reclutamiento = ceiling(total_dias),
    meses_reclutamiento = round(duracion_meses, 2)
  )
  
  return(resultados)
}
