# Patient-Specific Cardiovascular Digital Twins from ICU Signals

**Autora principal:** Diana Cabrera

**Coautor:** Luna, nombre completo pendiente de confirmación

**Institución:** Universidad Técnica del Norte, Ibarra, Ecuador

**Tipo de estudio:** computacional y metodológico

**Estado:** desarrollo metodológico; auditoría de datos completada y análisis exploratorio de ABP en progreso

## Propósito

Desarrollar y evaluar métodos reproducibles para caracterizar señales cardiovasculares de UCI y estudiar parámetros efectivos derivados de modelos fisiológicos.

Los objetivos actuales son:

1. caracterizar y controlar la calidad de señales de presión arterial invasiva (ABP);

2. detectar picos sistólicos y caracterizar latidos individuales;

3. estimar una constante de tiempo diastólica efectiva dependiente del modelo;

4. evaluar sensibilidad, incertidumbre e identificabilidad;

5. estudiar posteriormente asociaciones temporales con hipotensión o deterioro hemodinámico;

6. mantener separación por paciente y evitar fuga temporal durante la validación.

El proyecto no presenta el sistema como herramienta clínica ni presupone capacidad predictiva.

## Estado metodológico actual

El repositorio incluye:

\- lectura exploratoria de señales en formato WFDB;

\- control preliminar de calidad de ventanas ABP;

\- detección de picos sistólicos;

\- caracterización latido por latido;

\- estimación exploratoria de una constante de tiempo diastólica efectiva;

\- análisis de sensibilidad a la presión venosa asumida y al inicio del ajuste;

\- un modelo Windkessel de dos elementos para experimentos sintéticos con flujo conocido;

\- ocho pruebas automatizadas para los modelos sintéticos y las funciones ABP.

La constante de tiempo se interpreta como un parámetro efectivo dependiente del modelo. No se afirma que la resistencia vascular y la compliancia arterial sean identificables por separado utilizando únicamente ABP.

## Inicio rápido en Windows con Conda

```bat

conda env create -f environment.yml

conda activate cardiovascular-digital-twin

python -m pip install -e .

python -m pytest -q

