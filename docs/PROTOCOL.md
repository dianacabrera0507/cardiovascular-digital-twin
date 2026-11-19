# Protocolo maestro de investigación, versión 1.0

## 1. Información general
- Título provisional: *Patient-Specific Cardiovascular Digital Twins from ICU Signals: Physiological Parameter Estimation, Dynamic Modeling, and Early Hemodynamic Deterioration Analysis*.
- Autores: Diana Cabrera y Luna.
- Diseño: estudio computacional y metodológico.
- Estado: borrador; requiere completar campos tras revisión bibliográfica y auditoría.

## 2. Justificación y vacío de conocimiento
[PENDIENTE DE REVISIÓN BIBLIOGRÁFICA]

## 3. Objetivo general
Determinar si un modelo cardiovascular dinámico individualizado puede estimar parámetros fisiológicos interpretables a partir de señales de UCI y si su evolución temporal aporta información adicional para estudiar episodios de hipotensión o deterioro hemodinámico.

## 4. Objetivos específicos
1. Auditar disponibilidad, sincronización y calidad de señales.
2. Implementar un modelo Windkessel mínimo y extensiones justificadas.
3. Evaluar identificabilidad estructural y práctica.
4. Comparar modelos poblacionales e individualizados.
5. Evaluar estabilidad temporal e incertidumbre.
6. Estudiar asociaciones preevento sin fuga temporal.
7. Medir valor incremental frente a comparadores convencionales.
8. Cuantificar robustez y costo computacional.

## 5. Preguntas e hipótesis
Conservar las preguntas originales. Para cada una, especificar hipótesis nula, alternativa, estimando, población y criterio de decisión.

## 6. Fuente de datos
- MIMIC-IV: pendiente de confirmar versión, acceso, tablas y cobertura al ejecutar la auditoría.
- MIMIC-IV Waveform Database: pendiente de confirmar versión, registros y señales.
- Alternativa pública: solo tras verificar disponibilidad, licencia y adecuación.

## 7. Población, cohorte y unidad de análisis
[PENDIENTE DE AUDITORÍA DE DATOS]

### Inclusión
- Acceso autorizado.
- Enlace reproducible entre registro y unidad clínica, si aplica.
- Señales requeridas según el análisis.
- Duración y calidad mínimas preespecificadas.

### Exclusión
- Identificadores o sincronización insuficientes.
- Segmentos no interpretables según reglas predefinidas.
- Falta de información para el desenlace analizado.

## 8. Desenlaces
La definición primaria de hipotensión, umbral, duración, separación de episodios y análisis de sensibilidad se fijarán antes del análisis de resultados y se justificarán con literatura verificable.

## 9. Ventanas temporales
- Ventana de observación: [PENDIENTE]
- Gap de exclusión: [PENDIENTE]
- Horizonte: [PENDIENTE]
- Instante de predicción: explícito para cada muestra.

## 10. Señales y variables candidatas
ABP invasiva, NBP intermitente, ECG, PPG, HR, RR, SpO2 y variables clínicas adicionales. Su inclusión final depende de disponibilidad real. No tratar ondas y mediciones intermitentes como equivalentes.

## 11. Procesamiento
1. Inventario.
2. Verificación de unidades y frecuencia.
3. Alineación temporal.
4. Calidad y artefactos.
5. Filtrado documentado.
6. Detección de interrupciones.
7. Ventanas.
8. Registro de procedencia y transformaciones.

## 12. Modelo cardiovascular
Modelo Windkessel de partida:

C dP/dt = Qin - (P - Pv)/R

Las aproximaciones de Qin y Pv deberán justificarse. No separar R y C si los datos solo identifican una combinación como RC.

## 13. Estimación e identificabilidad
- Restricciones físicas.
- Múltiples inicializaciones.
- Perfiles de verosimilitud o costo.
- Sensibilidad local/global.
- Bootstrap o aproximación de incertidumbre.
- Recuperación en simulación.
- Validación sobre segmentos no usados en el ajuste.

## 14. Comparadores
- Windkessel poblacional.
- Windkessel individualizado.
- Parámetros por ventanas.
- Regresión/modelo clínico convencional.
- Árboles o boosting, si procede.
- Modelo temporal neuronal, si tamaño y diseño lo permiten.
- Híbrido fisiológico-neuronal, si demuestra viabilidad.

## 15. Métricas
### Ajuste
MAE, RMSE, error normalizado, residuos y estabilidad.

### Parámetros
Incertidumbre, variabilidad, sensibilidad, identificabilidad y recuperación sintética.

### Eventos
AUROC, AUPRC, sensibilidad, especificidad, calibración, falsas alarmas y anticipación.

### Eficiencia
Tiempo, memoria y número de parámetros.

## 16. Prevención de sesgos
- Separación por paciente.
- Ajuste de preprocesamiento dentro de cada partición.
- Variables disponibles antes de la predicción.
- Comparadores con información equivalente.
- Conjunto de prueba bloqueado.

## 17. Análisis de sensibilidad
Ruido, artefactos, faltantes, frecuencia, ventanas, filtros, definiciones de evento, regularización y canales disponibles.

## 18. Ética y gobernanza
No eludir acceso, no redistribuir datos restringidos y no enviar datos a servicios no autorizados.

## 19. Criterios de éxito y fracaso
Definir por objetivo. Un resultado negativo de identificabilidad o valor incremental se reportará sin reinterpretarlo como éxito predictivo.

## 20. Desviaciones del protocolo
Toda desviación se registrará con fecha, motivo, responsable e impacto esperado.
