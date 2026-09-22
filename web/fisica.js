/* Nucleo de calculo del receptor solar PPHE, portado a JavaScript.
 *
 * Traduccion literal de modelo_receptor/{mapa,props_fluido,receptor}.py y de
 * modelo_termohidraulico/{formulas,asignacion,geometria}.py. Mismas ecuaciones,
 * mismas constantes y misma biseccion anidada: lo que cambia es el lenguaje,
 * para que la pagina recalcule en el navegador sin servidor de Python.
 *
 * Todo en SI y las temperaturas en KELVIN, como en el original.
 */
(function (global) {
  'use strict';

  // ===========================================================================
  // Funcion error, que es la que sostiene todas las integrales del mapa
  // ===========================================================================

  // Aproximacion de Chebyshev de erfc (Numerical Recipes), precision ~1e-15
  // relativa. Hace falta porque JavaScript no trae erf en su libreria estandar.
  var ERFC_COF = [
    -1.3026537197817094, 6.4196979235649026e-1, 1.9476473204185836e-2,
    -9.561514786808631e-3, -9.46595344482036e-4, 3.66839497852761e-4,
    4.2523324806907e-5, -2.0278578112534e-5, -1.624290004647e-6,
    1.303655835580e-6, 1.5626441722e-8, -8.5238095915e-8, 6.529054439e-9,
    5.059343495e-9, -9.91364156e-10, -2.27365122e-10, 9.6467911e-11,
    2.394038e-12, -6.886027e-12, 8.94487e-13, 3.13092e-13, -1.12708e-13,
    3.81e-16, 7.106e-15
  ];

  function erfc(x) {
    var z = Math.abs(x);
    var t = 2.0 / (2.0 + z);
    var ty = 4.0 * t - 2.0;
    var d = 0.0, dd = 0.0, tmp;
    for (var j = ERFC_COF.length - 1; j > 0; j--) {
      tmp = d;
      d = ty * d - dd + ERFC_COF[j];
      dd = tmp;
    }
    var ans = t * Math.exp(-z * z + 0.5 * (ERFC_COF[0] + ty * d) - dd);
    return x >= 0.0 ? ans : 2.0 - ans;
  }

  function erf(x) { return 1.0 - erfc(x); }

  // ===========================================================================
  // Errores del modelo
  // ===========================================================================

  // Un fallo fisico (fluido fuera de rango, geometria sin correlacion) no es un
  // bug: es un punto de trabajo que no existe, y la pagina lo cuenta tal cual.
  function ErrorModelo(mensaje, consejo) {
    this.name = 'ErrorModelo';
    this.message = mensaje;
    this.consejo = consejo || '';
  }
  ErrorModelo.prototype = Object.create(Error.prototype);
  ErrorModelo.prototype.constructor = ErrorModelo;

  // ===========================================================================
  // Constantes del entorno y del material (receptor.py)
  // ===========================================================================

  var CERO_CELSIUS = 273.15;
  var R_UNIVERSAL = 8.314462618;
  var T_REF_ENTALPIA = 273.15;

  var SIGMA = 5.670374419e-8;   // Stefan-Boltzmann [W/(m2*K4)]
  var T_AMB = 298.15;           // Temperatura ambiente [K]
  var T_ENTORNO = 303.15;       // Entorno proximo, 30 C [K]
  var PESO_ENTORNO = 0.895;
  var PESO_CIELO = 0.955;

  // Ecuacion (6) del paper de M.R. Rodriguez
  var T_CIELO = Math.pow(
    (PESO_ENTORNO * Math.pow(T_ENTORNO, 4) + PESO_CIELO * Math.pow(T_AMB, 4)) /
    (PESO_ENTORNO + PESO_CIELO), 0.25);

  var LIMITE_AISI321 = 1150.0;  // Servicio continuo del AISI 321 [K]

  // Conductividad del AISI 321 [W/(m*K)]. Constante a proposito: la chapa es el
  // 2.7 % de la resistencia al fluido, y 20 corresponde a unos 700 K de chapa,
  // el centro del rango de trabajo. Ver la nota de receptor.py.
  var K_ACERO = 20.0;

  function perdidas(T_w, A, h_ext, eps) {
    return (h_ext * (T_w - T_AMB) +
            eps * SIGMA * (Math.pow(T_w, 4) - Math.pow(T_CIELO, 4))) * A;
  }

  // ===========================================================================
  // Propiedades del aire (props_fluido.py)
  // ===========================================================================

  function horner(coefs, x) {
    var r = 0.0;
    for (var i = coefs.length - 1; i >= 0; i--) r = r * x + coefs[i];
    return r;
  }

  var AIRE = {
    nombre: 'Aire',
    M: 28.96546e-3,
    coef_cp: [1086.4453, -659.81616, 1702.3009, -1486.6131, 586.96592, -88.252623],
    coef_k: [-0.25684666, 107.47369, -77.64382, 58.104571, -24.214681, 4.2141488],
    coef_mu: [0.45586184, 75.791245, -64.487066, 48.07116, -20.037935, 3.4861614],
    T_min: 250.0,
    T_max: 1600.0
  };
  AIRE.R = R_UNIVERSAL / AIRE.M;

  function verificarT(f, T) {
    if (T >= f.T_min && T <= f.T_max) return;
    throw new ErrorModelo(
      'El aire llega a ' + (T - CERO_CELSIUS).toFixed(0) + ' \u00b0C, fuera del ' +
      'rango de ajuste de las propiedades (' + (f.T_min - CERO_CELSIUS).toFixed(0) +
      ' a ' + (f.T_max - CERO_CELSIUS).toFixed(0) + ' \u00b0C).',
      'Sube el gasto G o baja el flujo de pico.');
  }

  function propiedades(f, T, p) {
    verificarT(f, T);
    var x = T / 1000.0;
    var rho = p / (f.R * T);
    var cp = horner(f.coef_cp, x);
    var k = horner(f.coef_k, x) * 1e-3;
    var mu = horner(f.coef_mu, x) * 1e-6;
    return { T: T, p: p, rho: rho, cp: cp, k: k, mu: mu, Pr: cp * mu / k };
  }

  function cpDe(f, T) { verificarT(f, T); return horner(f.coef_cp, T / 1000.0); }

  // Primitiva del polinomio de cp respecto de T, evaluada en x = T/1000.
  function integralCp(f, x) {
    var s = 0.0;
    for (var i = 0; i < f.coef_cp.length; i++) {
      s += f.coef_cp[i] * Math.pow(x, i + 1) / (i + 1);
    }
    return 1000.0 * s;
  }

  // Entalpia especifica referida a 273.15 K [J/kg]. Integral analitica del
  // polinomio de cp, no cp*(T-T_ref): con saltos de 200-400 K por columna, el
  // cp de la entrada se dejaria un 7 % del calor por el camino.
  function entalpia(f, T) {
    verificarT(f, T);
    return integralCp(f, T / 1000.0) - integralCp(f, T_REF_ENTALPIA / 1000.0);
  }

  // ===========================================================================
  // Geometrias de panel (geometria.py, Tabla 1 de O. Arsenyeva et al.)
  // ===========================================================================

  var PANELES = {
    PPHE1: { id: 'PPHE1', delta_pp: 0.8e-3, b_i: 3.4e-3, b: 5.5e-3, s_2l: 42e-3,
             s_t: 72e-3, d_sp: 7.2e-3, w_pp: 300e-3, l_pp: 1.0, w_e: 15e-3, Fx: 1.010 },
    PPHE2: { id: 'PPHE2', delta_pp: 1.0e-3, b_i: 3.0e-3, b: 1.5e-3, s_2l: 72e-3,
             s_t: 42e-3, d_sp: 7.2e-3, w_pp: 300e-3, l_pp: 1.0, w_e: 15e-3, Fx: 1.007 },
    PPHE3: { id: 'PPHE3', delta_pp: 1.0e-3, b_i: 7.0e-3, b: 20e-3, s_2l: 72e-3,
             s_t: 42e-3, d_sp: 7.2e-3, w_pp: 300e-3, l_pp: 1.0, w_e: 15e-3, Fx: 1.045 }
  };

  // ===========================================================================
  // Constantes de las leyes de potencia (asignacion.py, tablas de M. Piper)
  // ===========================================================================

  // Las tres familias de celda y sus rangos de validez, tal y como los dan las
  // tablas 2 y 3 de M. Piper: a = s2L/sT, b = dsp/sT, c = b_i/sT.
  var FAMILIAS = [
    { id: 'longitudinal', nombre: 'celda longitudinal',
      a: [0.57, 0.59], b: [0.10, 0.14], c: [0.042, 0.083],
      constantes: function (b, c) {
        return { n1: 8.74 * b + (17 * c + 0.73), n2: -0.38,
                 n3: 0.0775 * b + (0.38 * c + 0.005), n4: 0.75, n5: 0.4 };
      } },
    { id: 'cuadrada', nombre: 'celda cuadrada',
      a: [0.99, 1.01], b: [0.17, 0.24], c: [0.071, 0.143],
      constantes: function (b, c) {
        return { n1: -15.3 * b + (1.4 * c + 5.4), n2: 1.725 * b + (1.11 * c - 0.66),
                 n3: 0.03 * b + (0.76 * c - 0.032), n4: -1.12 * c + 0.905, n5: 0.4 };
      } },
    { id: 'transversal', nombre: 'celda transversal',
      a: [1.70, 1.72], b: [0.17, 0.24], c: [0.071, 0.170],
      constantes: function (b, c) {
        return { n1: 1.35 * b + (2.8 * c + 0.92), n2: 0.3 * b + (0.53 * c - 0.29),
                 n3: -0.163 * b + (0.711 * c + 0.022), n4: 0.29 * b + (-c + 0.8), n5: 0.4 };
      } }
  ];

  // Tolerancia relativa en los bordes de los rangos. Los pasos que elige la
  // pagina caen a veces justo en un borde (s_T = 72 mm da dsp/sT = 0.10), y el
  // redondeo de la division no puede dejarlos fuera por 1e-17.
  function dentro(v, r) { return v >= r[0] * (1 - 1e-9) && v <= r[1] * (1 + 1e-9); }

  function asignacionDeConstantes(sT, s2L, dsp, h) {
    var a = s2L / sT, b = dsp / sT, c = h / sT;
    for (var k = 0; k < FAMILIAS.length; k++) {
      var f = FAMILIAS[k];
      if (dentro(a, f.a) && dentro(b, f.b) && dentro(c, f.c)) {
        var n = f.constantes(b, c);
        n.familia = f.nombre;
        n.id = f.id;
        return n;
      }
    }
    throw new ErrorModelo(
      'La geometria del panel queda fuera de las correlaciones publicadas ' +
      '(s2L/sT = ' + a.toFixed(3) + ', dsp/sT = ' + b.toFixed(3) +
      ', b_i/sT = ' + c.toFixed(3) + ').',
      'Elige otro panel de la Tabla 1.');
  }

  // Lo que la correlacion deja elegir con un panel y una familia dados. d_sp y
  // b_i los fija el panel, asi que dsp/sT y b_i/sT acotan s_T, y s2L/sT acota
  // luego s_L. Devuelve los extremos en metros, o null si ningun s_T cumple a la
  // vez los dos rangos (PPHE3, con b_i = 7 mm, solo admite la transversal).
  function rangoGeometria(panel, idFamilia) {
    var f = null;
    for (var k = 0; k < FAMILIAS.length; k++) if (FAMILIAS[k].id === idFamilia) f = FAMILIAS[k];
    if (!f) return null;
    var sT_min = Math.max(panel.d_sp / f.b[1], panel.b_i / f.c[1]);
    var sT_max = Math.min(panel.d_sp / f.b[0], panel.b_i / f.c[0]);
    if (sT_min > sT_max * (1 + 1e-9)) return null;
    return { sT: [sT_min, sT_max], razon: f.a, familia: f };
  }

  // ===========================================================================
  // Formulas del canal interno (formulas.py)
  // ===========================================================================

  function deI(b_i) { return 2.0 * (b_i / Math.SQRT2); }           // Ec. (14)
  function reynolds(rho, u, dh, mu) { return rho * u * dh / mu; }
  function nusseltI(n3, n4, n5, Re, Pr) {                          // Ec. (20)
    return n3 * Math.pow(Re, n4) * Math.pow(Pr, n5);
  }
  function friccionI(n1, Re, n2) { return n1 * Math.pow(Re, n2); }  // Ec. (18)
  function fChI(b_i, w_pp, w_e) { return (b_i / Math.SQRT2) * (w_pp - 2 * w_e); }  // Ec. (16)

  // Resistencia hidraulica local de las zonas de distribucion, ec. (19). En el
  // receptor se reparte a medias entre la de entrada y la de salida.
  var ZETA_DZ = 1.5;

  // ===========================================================================
  // Mapa gaussiano de flujo incidente (mapa.py)
  // ===========================================================================

  var FWHM_POR_SIGMA = 2.0 * Math.sqrt(2.0 * Math.log(2.0));

  // Fraccion de campana normalizada que cae entre a y b. Escrita con erfc en
  // vez de con erf para no restar dos numeros casi iguales en las colas.
  function fraccionGaussiana(a, b, centro, sigma) {
    var t = sigma * Math.SQRT2;
    return 0.5 * (erfc((a - centro) / t) - erfc((b - centro) / t));
  }

  function Mapa(L, W, q_pico, sigma_x, sigma_y) {
    this.L = L; this.W = W;
    this.q_pico = q_pico;
    this.sigma_x = sigma_x; this.sigma_y = sigma_y;
  }

  Mapa.prototype.centro = function () { return [this.W / 2.0, this.L / 2.0]; };
  Mapa.prototype.area = function () { return this.L * this.W; };

  Mapa.prototype.q = function (x, y) {
    var xc = this.W / 2.0, yc = this.L / 2.0;
    return this.q_pico * Math.exp(
      -(Math.pow(x - xc, 2) / (2.0 * this.sigma_x * this.sigma_x) +
        Math.pow(y - yc, 2) / (2.0 * this.sigma_y * this.sigma_y)));
  };

  // Integral EXACTA sobre la placa, no sobre el plano infinito: la placa solo
  // recoge la parte central de la campana y el resto se derrama por los bordes.
  Mapa.prototype.potenciaTotal = function () {
    var xc = this.W / 2.0, yc = this.L / 2.0;
    return this.q_pico * (2.0 * Math.PI * this.sigma_x * this.sigma_y) *
      fraccionGaussiana(0.0, this.W, xc, this.sigma_x) *
      fraccionGaussiana(0.0, this.L, yc, this.sigma_y);
  };

  Mapa.prototype.qMedio = function () { return this.potenciaTotal() / this.area(); };

  Mapa.prototype.factorInterceptacion = function () {
    var xc = this.W / 2.0, yc = this.L / 2.0;
    return fraccionGaussiana(0.0, this.W, xc, this.sigma_x) *
           fraccionGaussiana(0.0, this.L, yc, this.sigma_y);
  };

  Mapa.prototype.factorPico = function () { return this.q_pico / this.qMedio(); };

  // Flujo MEDIO de cada celda: se INTEGRA la campana dentro de cada una, no se
  // muestrea en su centro, y asi la suma de potencias nodales reproduce
  // potenciaTotal() exactamente por bastas que sean las mallas.
  Mapa.prototype.mapaNodal = function (nx, ny) {
    var xc = this.W / 2.0, yc = this.L / 2.0;
    var dx = this.W / nx, dy = this.L / ny;
    var areaCelda = dx * dy;
    var fx = [], fy = [], i, j;
    for (i = 0; i < nx; i++) fx.push(fraccionGaussiana(i * dx, (i + 1) * dx, xc, this.sigma_x));
    for (j = 0; j < ny; j++) fy.push(fraccionGaussiana(j * dy, (j + 1) * dy, yc, this.sigma_y));
    var potenciaCampana = this.q_pico * 2.0 * Math.PI * this.sigma_x * this.sigma_y;
    var salida = [];
    for (j = 0; j < ny; j++) {
      var fila = [];
      for (i = 0; i < nx; i++) fila.push(potenciaCampana * fx[i] * fy[j] / areaCelda);
      salida.push(fila);
    }
    return salida;
  };

  // Despeja q_pico de la integral SOBRE LA PLACA, no de la del plano infinito.
  function mapaDesdePotencia(L, W, Q, sigma_x, sigma_y) {
    var provisional = new Mapa(L, W, 1.0, sigma_x, sigma_y);
    return new Mapa(L, W, Q / provisional.potenciaTotal(), sigma_x, sigma_y);
  }

  // ===========================================================================
  // Biseccion (receptor.py)
  // ===========================================================================

  function biseccion(f, a, b, tol, iteraciones) {
    tol = tol === undefined ? 1e-6 : tol;
    iteraciones = iteraciones === undefined ? 80 : iteraciones;
    var fa = f(a), fb = f(b);
    if (fa * fb > 0) {
      throw new ErrorModelo(
        'La biseccion no encierra la raiz en [' + a.toFixed(2) + ', ' + b.toFixed(2) + '].',
        'Punto de trabajo fuera del alcance del modelo.');
    }
    for (var n = 0; n < iteraciones; n++) {
      if (b - a < tol) break;
      var m = 0.5 * (a + b);
      var fm = f(m);
      if (fa * fm <= 0) { b = m; } else { a = m; fa = fm; }
    }
    return 0.5 * (a + b);
  }

  // ===========================================================================
  // El receptor
  // ===========================================================================

  /* cfg:
   *   L, W        [m]     alto y ancho de la placa
   *   q_pico      [W/m2]  flujo de pico  (o Q [W] si modo === 'potencia')
   *   Q           [W]
   *   modo        'pico' | 'potencia'
   *   sigma_x, sigma_y    [m]
   *   G           [kg/s]  gasto masico total
   *   T_ent       [K]     temperatura de entrada, uniforme
   *   p           [Pa]    presion a la entrada
   *   lam         [-]     mezcla entre columnas: 1 adiabaticas, 0 mezcla completa
   *   panel       id de PANELES
   *   s_t, s_2l   [m]     pasos de soldadura, opcionales: sustituyen a los del panel
   *   h_ext, eps, alfa
   *
   * La malla no se elige: la dicta el patron de soldaduras del panel, una
   * porcion por celda s_T x s_L. M = W/s_T y N = L/s_L, redondeados.
   */
  function resolver(cfg) {
    var panel = Object.assign({}, PANELES[cfg.panel] || PANELES.PPHE1);
    if (cfg.s_t) panel.s_t = cfg.s_t;
    if (cfg.s_2l) panel.s_2l = cfg.s_2l;
    var fluido = AIRE;

    var mapa = cfg.modo === 'potencia'
      ? mapaDesdePotencia(cfg.L, cfg.W, cfg.Q, cfg.sigma_x, cfg.sigma_y)
      : new Mapa(cfg.L, cfg.W, cfg.q_pico, cfg.sigma_x, cfg.sigma_y);

    var lam = cfg.lam;
    if (!(lam >= 0.0 && lam <= 1.0)) {
      throw new ErrorModelo('lambda tiene que estar entre 0 y 1.', '');
    }

    // -- Malla dictada por el patron de soldaduras ---------------------------
    var s_L = panel.s_2l / 2;
    var M = Math.max(Math.round(mapa.W / panel.s_t), 1);
    var N = Math.max(Math.round(mapa.L / s_L), 1);
    var dx = mapa.W / M, dy = mapa.L / N;
    if (dx <= panel.w_e || (M === 1 && mapa.W <= 2 * panel.w_e)) {
      throw new ErrorModelo(
        'Las soldaduras de borde (' + (panel.w_e * 1e3).toFixed(0) + ' mm) no caben ' +
        'en una columna de ' + (dx * 1e3).toFixed(0) + ' mm.',
        'Ensancha la placa.');
    }
    // Las columnas de los extremos pierden la soldadura de borde
    var anchos = [], i, j;
    for (i = 0; i < M; i++) anchos.push(dx);
    anchos[0] -= panel.w_e;
    anchos[M - 1] -= panel.w_e;
    var sumaAnchos = 0.0;
    for (i = 0; i < M; i++) sumaAnchos += anchos[i];
    // Un punto de soldadura por celda s_T x s_L (tresbolillo)
    var fracSoldadura = (Math.PI * panel.d_sp * panel.d_sp / 4) / (panel.s_t * s_L);

    var n = asignacionDeConstantes(panel.s_t, panel.s_2l, panel.d_sp, panel.b_i);
    var seccion = fChI(panel.b_i, mapa.W, panel.w_e);  // Ec. (16), w_pp = W
    var de = deI(panel.b_i);                            // Ec. (14)
    var e = panel.delta_pp;
    var h_ext = cfg.h_ext, eps = cfg.eps, alfa = cfg.alfa;
    var G = cfg.G, p = cfg.p, T_ent0 = cfg.T_ent;

    verificarT(fluido, T_ent0);

    // Coeficiente de pelicula del canal interno: ec. (3) -> Re, Pr -> Nu (20) -> h.
    // La velocidad es la misma en todas las columnas, porque el gasto se
    // reparte en proporcion a la anchura de paso.
    function hInterno(T) {
      var pr = propiedades(fluido, T, p);
      var u = G / (pr.rho * seccion);
      var Re = reynolds(pr.rho, u, de, pr.mu);
      var Nu = nusseltI(n.n3, n.n4, n.n5, Re, pr.Pr);
      return { h: Nu * pr.k / de, u: u, Re: Re, Pr: pr.Pr, Nu: Nu };
    }

    // Perdida de carga por friccion por unidad de longitud, a T y p [Pa/m].
    // Primer termino de la ec. (19), con f de la ec. (18).
    function friccion(T, pLocal) {
      var pr = propiedades(fluido, T, pLocal);
      var u = G / (pr.rho * seccion);
      var Re = reynolds(pr.rho, u, de, pr.mu);
      return friccionI(n.n1, Re, n.n2) / de * pr.rho * u * u / 2;
    }

    var A = dx * dy;
    var G_col = [], A_f = [];
    for (i = 0; i < M; i++) {
      G_col.push(G * anchos[i] / sumaAnchos);
      A_f.push(anchos[i] * dy * (1.0 - fracSoldadura));
    }

    // -- Balance de una porcion ---------------------------------------------
    // A recibe el flujo y pierde calor; A_f, sin soldaduras, lo cede al fluido.
    function porcion(q_inc, T_ent, Gc, Af) {
      var q_abs = alfa * q_inc;

      function cerrar(T_sal) {
        var T_f = 0.5 * (T_ent + T_sal);
        var h_int = hInterno(T_f).h;

        // Con K_ACERO constante la conductancia ya no depende de T_w: se
        // calcula aqui una vez, y no en cada evaluacion de la biseccion.
        var UA = Af / (e / K_ACERO + 1.0 / h_int);

        function balancePared(T_w) {
          return q_abs - perdidas(T_w, A, h_ext, eps) - UA * (T_w - T_f);
        }

        // La pared esta entre el cielo (pierde mas de lo que recibe) y la
        // temperatura a la que solo la radiacion ya se lleva todo el flujo.
        var T_w_max = Math.max(
          T_f, Math.pow(q_abs / (A * eps * SIGMA) + Math.pow(T_CIELO, 4), 0.25)) + 1.0;
        var T_w = biseccion(balancePared, T_CIELO, T_w_max);
        return { q: UA * (T_w - T_f), T_w: T_w };
      }

      function desequilibrio(T_sal) {
        return cerrar(T_sal).q - Gc * (entalpia(fluido, T_sal) - entalpia(fluido, T_ent));
      }

      // Cota superior: todo el calor absorbido al fluido con el cp de la
      // entrada, que es el menor del tramo. Como el cp del aire crece con T, la
      // temperatura real queda por debajo.
      var T_max = Math.min(T_ent + q_abs / (Gc * cpDe(fluido, T_ent)), fluido.T_max);
      if (desequilibrio(T_max) > 0) {
        throw new ErrorModelo(
          'El aire se saldria de ' + (fluido.T_max - CERO_CELSIUS).toFixed(0) +
          ' °C dentro de una porcion.',
          'Sube el gasto G o baja el flujo incidente.');
      }
      var T_sal = biseccion(desequilibrio, fluido.T_min, T_max);
      var fin = cerrar(T_sal);
      return { T_sal: T_sal, T_w: fin.T_w, q: fin.q };
    }

    // -- Mezcla entre columnas ----------------------------------------------
    // Pondera en ENTALPIA el caso adiabatico y el de mezcla completa, para que
    // la mezcla conserve la energia para cualquier lambda.
    function mezclar(T_ad) {
      var T_bajo = Math.min.apply(null, T_ad), T_alto = Math.max.apply(null, T_ad);
      if (lam === 1.0 || T_alto - T_bajo < 1e-9) return T_ad.slice();
      var h_ad = [], h_mez = 0.0, k;
      for (k = 0; k < M; k++) {
        h_ad.push(entalpia(fluido, T_ad[k]));
        h_mez += G_col[k] * h_ad[k];
      }
      h_mez /= G;
      var a = Math.max(T_bajo - 1.0, fluido.T_min), b = Math.min(T_alto + 1.0, fluido.T_max);
      var salida = [];
      for (k = 0; k < M; k++) {
        var objetivo = lam * h_ad[k] + (1.0 - lam) * h_mez;
        salida.push(biseccion(function (T) { return entalpia(fluido, T) - objetivo; },
                              a, b, 1e-9));
      }
      return salida;
    }

    // -- Recorrido de la placa, fila a fila y de abajo arriba ----------------
    var flujo = mapa.mapaNodal(M, N);
    var T_fluido = [], T_adiab = [], T_pared = [], q_nodo = [];
    var q_util = 0.0;
    var T = [];
    for (i = 0; i < M; i++) T.push(T_ent0);
    for (j = 0; j < N; j++) {
      var T_ad = new Array(M), T_w = new Array(M), qf = new Array(M);
      for (i = 0; i < M; i++) {
        var r = porcion(flujo[j][i] * A, T[i], G_col[i], A_f[i]);
        T_ad[i] = r.T_sal; T_w[i] = r.T_w; qf[i] = r.q;
        q_util += r.q;
      }
      T = mezclar(T_ad);
      T_fluido.push(T); T_adiab.push(T_ad); T_pared.push(T_w); q_nodo.push(qf);
    }

    // -- Perdida de carga, columna a columna y bajando la presion ------------
    // El calculo termico no depende de la presion (Re = G''*de/mu), asi que se
    // hace despues y aparte. Cada densidad se evalua a la presion local.
    var G2 = Math.pow(G / seccion, 2);
    var dpFriccion = [], dpDistribucion = [], dpAceleracion = [];
    for (i = 0; i < M; i++) {
      var Tc = T_ent0, pc = p;
      var v = 1.0 / propiedades(fluido, Tc, pc).rho;
      var fr = 0.0, ac = 0.0;
      var di = ZETA_DZ / 2 * G2 * v;                 // Zona de distribucion de entrada
      pc -= di;
      for (j = 0; j < N; j++) {
        var T_s = T_fluido[j][i];
        var dpf = friccion(0.5 * (Tc + T_s), pc) * dy;
        var v_s = 1.0 / propiedades(fluido, T_s, pc - dpf).rho;
        var dpa = G2 * (v_s - v);
        fr += dpf; ac += dpa;
        pc -= dpf + dpa;
        if (pc <= 0.0) {
          throw new ErrorModelo(
            'La perdida de carga agota los ' + (p / 1e5).toFixed(1) +
            ' bar de entrada antes de la salida.',
            'Sube la presion o baja el gasto.');
        }
        Tc = T_s; v = v_s;
      }
      di += ZETA_DZ / 2 * G2 * v;                    // Zona de distribucion de salida
      dpFriccion.push(fr); dpDistribucion.push(di); dpAceleracion.push(ac);
    }
    function mediaGasto(valores) {
      var s = 0.0;
      for (var k = 0; k < M; k++) s += G_col[k] * valores[k];
      return s / G;
    }
    var dpColumnas = [];
    for (i = 0; i < M; i++) dpColumnas.push(dpFriccion[i] + dpDistribucion[i] + dpAceleracion[i]);

    // -- Resumen -------------------------------------------------------------
    var perfilSalida = T_fluido[N - 1];
    // Temperatura de mezcla a la salida, por entalpia
    var h_media = entalpia(fluido, T_ent0) + q_util / G;
    var T_media_sal = biseccion(function (Tx) { return entalpia(fluido, Tx) - h_media; },
                                fluido.T_min, fluido.T_max, 1e-9);

    var canal = hInterno(0.5 * (T_ent0 + T_media_sal));

    var T_pared_max = -Infinity, T_pared_min = Infinity;
    var T_fluido_max = -Infinity, T_fluido_min = Infinity;
    for (j = 0; j < N; j++) for (i = 0; i < M; i++) {
      if (T_pared[j][i] > T_pared_max) T_pared_max = T_pared[j][i];
      if (T_pared[j][i] < T_pared_min) T_pared_min = T_pared[j][i];
      if (T_fluido[j][i] > T_fluido_max) T_fluido_max = T_fluido[j][i];
      if (T_fluido[j][i] < T_fluido_min) T_fluido_min = T_fluido[j][i];
    }

    var Q = mapa.potenciaTotal();
    var xCentros = [], yCentros = [];
    for (i = 0; i < M; i++) xCentros.push(mapa.W * (i + 0.5) / M);
    for (j = 0; j < N; j++) yCentros.push(mapa.L * (j + 0.5) / N);

    return {
      mapa: mapa, panel: panel, n: n, M: M, N: N, dx: dx, dy: dy, lam: lam,
      flujo: flujo, T_fluido: T_fluido, T_adiabatica: T_adiab, T_pared: T_pared,
      q_nodo: q_nodo,
      A_celda: A, A_f: A_f, G_col: G_col, fracSoldadura: fracSoldadura, T_ent: T_ent0,
      perfilSalida: perfilSalida, xCentros: xCentros, yCentros: yCentros,
      Q_incidente: Q, q_util: q_util, perdidas: Q - q_util,
      rendimiento: q_util / Q,
      T_media_salida: T_media_sal,
      T_max_salida: Math.max.apply(null, perfilSalida),
      T_min_salida: Math.min.apply(null, perfilSalida),
      T_pared_max: T_pared_max, T_pared_min: T_pared_min,
      T_fluido_max: T_fluido_max, T_fluido_min: T_fluido_min,
      q_pico: mapa.q_pico, q_medio: mapa.qMedio(),
      interceptacion: mapa.factorInterceptacion(),
      factorPico: mapa.factorPico(),
      u: canal.u, Re: canal.Re, h_int: canal.h, Pr: canal.Pr, Nu: canal.Nu,
      de: de, seccion: seccion,
      dpColumnas: dpColumnas,
      perdidaCarga: mediaGasto(dpColumnas),
      dpFriccion: mediaGasto(dpFriccion),
      dpDistribucion: mediaGasto(dpDistribucion),
      dpAceleracion: mediaGasto(dpAceleracion),
      limiteSuperado: T_pared_max > LIMITE_AISI321
    };
  }

  // ===========================================================================

  var RCP = {
    erf: erf, erfc: erfc,
    CERO_CELSIUS: CERO_CELSIUS, SIGMA: SIGMA, T_AMB: T_AMB, T_CIELO: T_CIELO,
    LIMITE_AISI321: LIMITE_AISI321, FWHM_POR_SIGMA: FWHM_POR_SIGMA,
    AIRE: AIRE, PANELES: PANELES, ErrorModelo: ErrorModelo,
    propiedades: propiedades, entalpia: entalpia, K_ACERO: K_ACERO, ZETA_DZ: ZETA_DZ,
    LAMBDA: 0.7,
    asignacionDeConstantes: asignacionDeConstantes,
    FAMILIAS: FAMILIAS, rangoGeometria: rangoGeometria,
    Mapa: Mapa, mapaDesdePotencia: mapaDesdePotencia,
    resolver: resolver
  };

  global.RCP = RCP;
  if (typeof module !== 'undefined' && module.exports) module.exports = RCP;
})(typeof globalThis !== 'undefined' ? globalThis : this);
