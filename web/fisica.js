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

  function asignacionDeConstantes(sT, s2L, dsp, h) {
    var a = s2L / sT, b = dsp / sT, c = h / sT;
    if (a >= 0.57 && a <= 0.59 && b >= 0.10 && b <= 0.14 && c >= 0.042 && c <= 0.083) {
      return { familia: 'celda longitudinal',
               n1: 8.74 * b + (17 * c + 0.73), n2: -0.38,
               n3: 0.0775 * b + (0.38 * c + 0.005), n4: 0.75, n5: 0.4 };
    }
    if (a >= 0.99 && a <= 1.01 && b >= 0.17 && b <= 0.24 && c >= 0.071 && c <= 0.143) {
      return { familia: 'celda cuadrada',
               n1: -15.3 * b + (1.4 * c + 5.4), n2: 1.725 * b + (1.11 * c - 0.66),
               n3: 0.03 * b + (0.76 * c - 0.032), n4: -1.12 * c + 0.905, n5: 0.4 };
    }
    if (a >= 1.70 && a <= 1.72 && b >= 0.17 && b <= 0.24 && c >= 0.071 && c <= 0.170) {
      return { familia: 'celda transversal',
               n1: 1.35 * b + (2.8 * c + 0.92), n2: 0.3 * b + (0.53 * c - 0.29),
               n3: -0.163 * b + (0.711 * c + 0.022), n4: 0.29 * b + (-c + 0.8), n5: 0.4 };
    }
    throw new ErrorModelo(
      'La geometria del panel queda fuera de las correlaciones publicadas ' +
      '(s2L/sT = ' + a.toFixed(3) + ', dsp/sT = ' + b.toFixed(3) +
      ', b_i/sT = ' + c.toFixed(3) + ').',
      'Elige otro panel de la Tabla 1.');
  }

  // ===========================================================================
  // Formulas del canal interno (formulas.py)
  // ===========================================================================

  function deI(b_i) { return 2.0 * (b_i / Math.SQRT2); }           // Ec. (14)
  function reynolds(rho, u, dh, mu) { return rho * u * dh / mu; }
  function nusseltI(n3, n4, n5, Re, Pr) {                          // Ec. (20)
    return n3 * Math.pow(Re, n4) * Math.pow(Pr, n5);
  }

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
   *   p           [Pa]
   *   N, M                porciones en altura y en anchura
   *   panel       id de PANELES
   *   h_ext, eps, alfa
   */
  function resolver(cfg) {
    var panel = PANELES[cfg.panel] || PANELES.PPHE1;
    var fluido = AIRE;

    var mapa = cfg.modo === 'potencia'
      ? mapaDesdePotencia(cfg.L, cfg.W, cfg.Q, cfg.sigma_x, cfg.sigma_y)
      : new Mapa(cfg.L, cfg.W, cfg.q_pico, cfg.sigma_x, cfg.sigma_y);

    var n = asignacionDeConstantes(panel.s_t, panel.s_2l, panel.d_sp, panel.b_i);
    var seccion = panel.b_i / Math.SQRT2 * mapa.W;   // Ec. (16), sin soldaduras de borde
    var de = deI(panel.b_i);                          // Ec. (14)
    var e = panel.delta_pp;
    var h_ext = cfg.h_ext, eps = cfg.eps, alfa = cfg.alfa;
    var G = cfg.G, p = cfg.p, T_ent0 = cfg.T_ent;
    var N = cfg.N, M = cfg.M;

    verificarT(fluido, T_ent0);

    // Coeficiente de pelicula del canal interno: ec. (3) -> Re, Pr -> Nu (20) -> h.
    // No depende de M: el gasto de una columna y su seccion de paso son ambos
    // proporcionales a dx, y el cociente se cancela.
    function hInterno(T) {
      var pr = propiedades(fluido, T, p);
      var u = G / (pr.rho * seccion);
      var Re = reynolds(pr.rho, u, de, pr.mu);
      var Nu = nusseltI(n.n3, n.n4, n.n5, Re, pr.Pr);
      return { h: Nu * pr.k / de, u: u, Re: Re, Pr: pr.Pr, Nu: Nu };
    }

    var A = (mapa.W / M) * (mapa.L / N);
    var G_col = G / M;

    // -- Balance de una porcion ---------------------------------------------
    function porcion(q_inc, T_ent) {
      var q_abs = alfa * q_inc;

      function cerrar(T_sal) {
        var T_f = 0.5 * (T_ent + T_sal);
        var h_int = hInterno(T_f).h;

        // Con K_ACERO constante la resistencia ya no depende de T_w: se calcula
        // aqui una vez, y no en cada evaluacion de la biseccion de abajo.
        var resistencia = e / K_ACERO + 1.0 / h_int;

        function balancePared(T_w) {
          return q_abs - perdidas(T_w, A, h_ext, eps) - A * (T_w - T_f) / resistencia;
        }

        // La pared esta entre el cielo (pierde mas de lo que recibe) y la
        // temperatura a la que solo la radiacion ya se lleva todo el flujo.
        var T_w_max = Math.max(
          T_f, Math.pow(q_abs / (A * eps * SIGMA) + Math.pow(T_CIELO, 4), 0.25)) + 1.0;
        var T_w = biseccion(balancePared, T_CIELO, T_w_max);
        return { q: A * (T_w - T_f) / resistencia, T_w: T_w };
      }

      function desequilibrio(T_sal) {
        return cerrar(T_sal).q - G_col * (entalpia(fluido, T_sal) - entalpia(fluido, T_ent));
      }

      // Cota superior: todo el calor absorbido al fluido con el cp de la
      // entrada, que es el menor del tramo. Como el cp del aire crece con T, la
      // temperatura real queda por debajo.
      var T_max = Math.min(T_ent + q_abs / (G_col * cpDe(fluido, T_ent)), fluido.T_max);
      if (desequilibrio(T_max) > 0) {
        throw new ErrorModelo(
          'El aire se saldria de ' + (fluido.T_max - CERO_CELSIUS).toFixed(0) +
          ' \u00b0C dentro de una porcion.',
          'Sube el gasto G o baja el flujo incidente.');
      }
      var T_sal = biseccion(desequilibrio, fluido.T_min, T_max);
      var fin = cerrar(T_sal);
      return { T_sal: T_sal, T_w: fin.T_w, q: fin.q };
    }

    // -- Recorrido de la placa, columna a columna y de abajo arriba ----------
    var flujo = mapa.mapaNodal(M, N);
    var T_fluido = [], T_pared = [], q_nodo = [], j, i;
    for (j = 0; j < N; j++) {
      T_fluido.push(new Array(M)); T_pared.push(new Array(M)); q_nodo.push(new Array(M));
    }

    var q_util = 0.0;
    for (i = 0; i < M; i++) {
      var T = T_ent0;
      for (j = 0; j < N; j++) {
        var r = porcion(flujo[j][i] * A, T);
        T = r.T_sal;
        T_fluido[j][i] = T;
        T_pared[j][i] = r.T_w;
        q_nodo[j][i] = r.q;
        q_util += r.q;
      }
    }

    // -- Resumen -------------------------------------------------------------
    var perfilSalida = T_fluido[N - 1];
    var T_media_sal = 0.0;
    for (i = 0; i < M; i++) T_media_sal += perfilSalida[i];
    T_media_sal /= M;

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
      mapa: mapa, panel: panel, n: n,
      flujo: flujo, T_fluido: T_fluido, T_pared: T_pared, q_nodo: q_nodo,
      A_celda: A, G_col: G_col, T_ent: T_ent0,
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
      limiteSuperado: T_pared_max > LIMITE_AISI321
    };
  }

  // ===========================================================================

  var RCP = {
    erf: erf, erfc: erfc,
    CERO_CELSIUS: CERO_CELSIUS, SIGMA: SIGMA, T_AMB: T_AMB, T_CIELO: T_CIELO,
    LIMITE_AISI321: LIMITE_AISI321, FWHM_POR_SIGMA: FWHM_POR_SIGMA,
    AIRE: AIRE, PANELES: PANELES, ErrorModelo: ErrorModelo,
    propiedades: propiedades, entalpia: entalpia, K_ACERO: K_ACERO,
    asignacionDeConstantes: asignacionDeConstantes,
    Mapa: Mapa, mapaDesdePotencia: mapaDesdePotencia,
    resolver: resolver
  };

  global.RCP = RCP;
  if (typeof module !== 'undefined' && module.exports) module.exports = RCP;
})(typeof globalThis !== 'undefined' ? globalThis : this);
