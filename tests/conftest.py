"""
Configuracion comun de pytest.

Anade al ``sys.path`` la raiz del repositorio (para importar ``simulation``)
y ``src/`` (para importar ``circust``), de modo que la suite corra sin
instalar el paquete.
"""
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent

for extra in (PROJECT_ROOT, PROJECT_ROOT / "src"):
    if str(extra) not in sys.path:
        sys.path.insert(0, str(extra))
