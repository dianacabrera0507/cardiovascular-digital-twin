# Patient-Specific Cardiovascular Digital Twins from ICU Signals

**Autores:** Diana Cabrera y Luna  
**Tipo de estudio:** computacional y metodológico  
**Estado:** estructura inicial, sin datos clínicos ni resultados

## Propósito

Desarrollar y evaluar un modelo cardiovascular dinámico específico por paciente para:

1. representar señales cardiovasculares de UCI;
2. estimar parámetros latentes interpretables del modelo;
3. evaluar identificabilidad, incertidumbre y estabilidad;
4. estudiar asociaciones temporales con hipotensión o deterioro hemodinámico;
5. comparar modelos fisiológicos, estadísticos, neuronales e híbridos.

El proyecto no presenta el sistema como herramienta clínica ni presupone capacidad predictiva.

## Inicio rápido

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
pip install -e .
pytest
python scripts/run_synthetic_demo.py --config configs/base.yaml
```

## Datos

Los datos MIMIC/PhysioNet no se incluyen en este repositorio. Coloque exclusivamente datos obtenidos y procesados conforme a los acuerdos aplicables. Nunca confirme archivos clínicos restringidos en Git.

- `data/raw/`: solo lectura; datos originales autorizados.
- `data/interim/`: resultados temporales.
- `data/processed/`: cohortes y ventanas derivadas autorizadas.
- `data/synthetic/`: datos sintéticos para desarrollo y pruebas.

## Flujo científico

1. Auditoría de datos.
2. Cohorte y eventos preespecificados.
3. Control de calidad y sincronización.
4. Simulación y recuperación de parámetros.
5. Identificabilidad estructural y práctica.
6. Modelos poblacional, individualizado, dinámico e híbrido.
7. Evaluación fuera de muestra.
8. Asociación y predicción sin fuga temporal.
9. Robustez, incertidumbre, calibración y eficiencia.
10. Manuscrito y material suplementario.

## Reglas de interpretación

- Reconstrucción, estimación de parámetros y predicción son objetivos distintos.
- Los parámetros sin referencia independiente se denominan **parámetros estimados del modelo**.
- La partición principal debe realizarse por paciente.
- Toda variable predictora debe estar disponible antes del instante de predicción.
- No se inventan tamaños de cohorte, resultados ni variables disponibles.

Consulte `docs/PROTOCOL.md`, `docs/ANALYSIS_PLAN.md` y `docs/DATA_GOVERNANCE.md`.
