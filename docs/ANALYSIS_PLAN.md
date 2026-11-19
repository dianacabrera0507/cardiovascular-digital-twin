# Plan de análisis

## Fase A: auditoría
- Crear inventario por paciente, estancia, registro y canal.
- Resumir duración, frecuencia, unidades, interrupciones y faltantes.
- Generar diagrama de inclusión sin anticipar cifras.

## Fase B: simulación
- Simular parámetros conocidos.
- Añadir ruido, artefactos, faltantes y submuestreo.
- Medir sesgo, varianza, RMSE y cobertura.
- Examinar si R, C o RC son recuperables.

## Fase C: reconstrucción
- Ajustar en una ventana.
- Evaluar en una ventana temporal no usada.
- Comparar modelo común e individualizado.

## Fase D: análisis longitudinal
- Estimar nivel, pendiente, variabilidad, cambios e incertidumbre.
- Usar modelos compatibles con medidas repetidas por paciente.

## Fase E: eventos
- Bloquear definición, horizonte y gap antes de evaluar.
- Comparar variables convencionales, parámetros y combinación.
- Reportar discriminación, calibración y falsas alarmas.

## Fase F: robustez
- Ablaciones por señales y componentes.
- Sensibilidad a preprocesamiento y definición del evento.
- Evaluación del costo computacional.

## Multiplicidad
Predefinir análisis primario, secundarios y exploratorios. Aplicar corrección cuando corresponda e interpretar intervalos de confianza junto con tamaños de efecto.
