/* Nucleo de calculo del receptor octogonal de sal solar, portado a JavaScript.
 *
 * Traduccion de modelo_receptor/receptor_octogonal.py, con lo que usa de
 * receptor.py, props_fluido.py (SAL_SOLAR) y modelo_termohidraulico/{formulas,
 * asignacion}.py. Mismas ecuaciones, mismas constantes y el mismo resolvedor:
 * Newton con salvaguarda de intervalo en la salida de cada porcion, Newton
 * monotono en la pared, y punto fijo del gasto acelerado con la secante.
 *
 * La unica diferencia numerica es de detalle: numpy resolvia las M columnas de
 * una fila a la vez e iteraba hasta que convergia la ultima; aqui cada porcion
 * para cuando converge ella. Las dos cumplen la misma tolerancia.
 *
 * Incluye la estrella de 8 puntas (cfg.alfa_estrella): placas dobladas de
 * ancho w/cos(alfa) con flujo cos(alfa)*q'', como en el original.
 *
 * Todo en SI y las temperaturas en KELVIN, como en el original.
 */
(function (global) {
  'use strict';

  // ===========================================================================
  // Funcion error (Chebyshev de erfc, Numerical Recipes, ~1e-15 relativa)
  // ===========================================================================

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

  // ===========================================================================
  // Errores del modelo: un punto de trabajo que no existe, no un fallo
  // ===========================================================================

  function ErrorModelo(mensaje, consejo) {
    this.name = 'ErrorModelo';
    this.message = mensaje;
    this.consejo = consejo || '';
  }
  ErrorModelo.prototype = Object.create(Error.prototype);
  ErrorModelo.prototype.constructor = ErrorModelo;

  // ===========================================================================
  // Constantes (receptor.py y receptor_octogonal.py)
  // ===========================================================================

  var CERO_CELSIUS = 273.15;
  var SIGMA = 5.670374419e-8;
  var T_AMB = 298.15;
  var T_ENTORNO = 303.15;
  var PESO_ENTORNO = 0.895;
  var PESO_CIELO = 0.955;
  var T_CIELO = Math.pow(
    (PESO_ENTORNO * Math.pow(T_ENTORNO, 4) + PESO_CIELO * Math.pow(T_AMB, 4)) /
    (PESO_ENTORNO + PESO_CIELO), 0.25);

  var K_ACERO = 20.0;                    // AISI 321, constante [W/(m*K)]
  var LIMITE_AISI321 = 1150.0;           // Servicio continuo del acero [K]
  var T_PELICULA_MAX = CERO_CELSIUS + 600.0;  // Descomposicion de la sal [K]
  var GRAVEDAD = 9.80665;
  var ZETA_DZ = 1.5;                     // Zonas de distribucion, ec. (19)
  var TOL_FILA = 1e-9;                   // Newton en la salida de cada porcion [K]
  var N_PLACAS = 8;
  // Angulo de la estrella con el que los picos internos llegan al centro
  var ALFA_ESTRELLA_MAX = Math.PI / 2 - Math.PI / N_PLACAS;

  // Rangos de validez de las correlaciones, tal y como los dan los articulos
  // (Fuentes/ del TFG):
  //  - M. Piper et al., Int. J. Therm. Sci. 120 (2017) 459-468, tablas 2, 3 y
  //    4: friccion (ec. 2) y Nusselt (ec. 6) para 1000 <= Re <= 8000 y
  //    1 <= Pr <= 150. Las constantes n1..n5 del modelo son las de esas tablas,
  //    asi que este es el rango que gobierna el calculo.
  //  - O. Arsenyeva et al., Appl. Therm. Eng. 147 (2019) 579-591: las ecs. (20)
  //    y (21), y la friccion del canal externo, para 9500 <= Re <= 30000.
  var RANGOS = {
    piper: { nombre: 'Piper et al. (2017)', Re: [1000, 8000], Pr: [1, 150],
             ecuaciones: 'fricción y el Nusselt del canal interno (ecs. 2 y 6, tablas 2 y 3)' },
    arsenyeva: { nombre: 'Arsenyeva et al. (2019)', Re: [9500, 30000], Pr: null,
                 ecuaciones: 'ecs. (20) y (21)' }
  };

  // round() de Python: los empates van al PAR (262.5 -> 262), no hacia arriba
  // como Math.round. La malla es M = round(w/s_T) y N = round(L/s_2l), y con
  // pasos redondos (L = 10.5 m, s_2l = 40 mm) el empate exacto se da.
  function redondeoPython(x) {
    var r = Math.round(x);
    if (Math.abs(x - Math.trunc(x)) === 0.5 && r % 2 !== 0) r -= 1;
    return r;
  }

  function perdidas(T_w, A, h_ext, eps) {
    var T2 = T_w * T_w, Tc2 = T_CIELO * T_CIELO;
    return (h_ext * (T_w - T_AMB) + eps * SIGMA * (T2 * T2 - Tc2 * Tc2)) * A;
  }

  // ===========================================================================
  // Sal solar (props_fluido.py, correlaciones de Zavoico, t en C)
  // ===========================================================================

  function horner(coefs, x) {
    var r = 0.0;
    for (var i = coefs.length - 1; i >= 0; i--) r = r * x + coefs[i];
    return r;
  }

  var SAL = {
    nombre: 'Sal solar',
    coef_rho: [2090.0, -0.636],
    coef_cp: [1443.0, 0.172],
    coef_k: [0.443, 1.9e-4],
    coef_mu: [22.714e-3, -0.120e-3, 2.281e-7, -1.474e-10],
    T_min: CERO_CELSIUS + 260.0,
    T_max: CERO_CELSIUS + 600.0
  };

  function verificarT(T) {
    if (T >= SAL.T_min && T <= SAL.T_max) return;
    var fuera = T > SAL.T_max;
    throw new ErrorModelo(
      'La sal llega a ' + (T - CERO_CELSIUS).toFixed(0) + ' °C, fuera del rango de ' +
      'sus propiedades (' + (SAL.T_min - CERO_CELSIUS).toFixed(0) + ' a ' +
      (SAL.T_max - CERO_CELSIUS).toFixed(0) + ' °C).',
      fuera ? 'Baja la temperatura de salida o baja λ para que las columnas se mezclen más.'
            : 'Las pérdidas superan al flujo absorbido: sube el flujo o la absortividad.');
  }

  function propiedades(T) {
    verificarT(T);
    var t = T - CERO_CELSIUS;
    var cp = horner(SAL.coef_cp, t), k = horner(SAL.coef_k, t), mu = horner(SAL.coef_mu, t);
    return { T: T, rho: horner(SAL.coef_rho, t), cp: cp, k: k, mu: mu, Pr: cp * mu / k };
  }
  function rho(T) { verificarT(T); return horner(SAL.coef_rho, T - CERO_CELSIUS); }
  function cp(T) { verificarT(T); return horner(SAL.coef_cp, T - CERO_CELSIUS); }

  // Primitiva exacta del polinomio de cp, referida a 0 C
  function hSinVerificar(t) {
    var s = 0.0;
    for (var i = 0; i < SAL.coef_cp.length; i++) s += SAL.coef_cp[i] * Math.pow(t, i + 1) / (i + 1);
    return s;
  }
  function entalpia(T) { verificarT(T); return hSinVerificar(T - CERO_CELSIUS); }

  // Inversa de h(T) por Newton, desde la estimacion con cp constante
  function TdesdeH(h) {
    var T = CERO_CELSIUS + h / SAL.coef_cp[0];
    for (var it = 0; it < 30; it++) {
      var t = T - CERO_CELSIUS;
      var paso = (hSinVerificar(t) - h) / horner(SAL.coef_cp, t);
      T -= paso;
      if (Math.abs(paso) < 1e-10) break;
    }
    verificarT(T);
    return T;
  }

  // ===========================================================================
  // Constantes de Piper (asignacion.py) y formulas del canal (formulas.py)
  // ===========================================================================

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
      'La geometría del panel queda fuera de las correlaciones de Piper ' +
      '(s2L/sT = ' + a.toFixed(3) + ', dsp/sT = ' + b.toFixed(3) + ', b_i/sT = ' + c.toFixed(3) + ').',
      'Elige otra combinación de pasos.');
  }

  function deI(b_i) { return 2.0 * (b_i / Math.SQRT2); }                        // Ec. (14)
  function fChI(b_i, w_pp, w_e) { return (b_i / Math.SQRT2) * (w_pp - 2 * w_e); } // Ec. (16)

  // ===========================================================================
  // Mapa de la mitad del receptor (MapaMitad)
  // ===========================================================================

  function fraccionGaussiana(a, b, centro, sigma) {
    var t = sigma * Math.SQRT2;
    return 0.5 * (erfc((a - centro) / t) - erfc((b - centro) / t));
  }

  // Gaussiana del receptor desarrollado vista sobre la mitad x en [0, W]: pico
  // en el plano de simetria (x = 0) y a media altura.
  function MapaMitad(L, w, n, q_pico, sx, sy) {
    this.L = L; this.w = w; this.n = n; this.q_pico = q_pico;
    this.sigma_x = sx; this.sigma_y = sy; this.W = n * w;
  }
  MapaMitad.prototype.fx = function (a, b) { return fraccionGaussiana(a, b, 0.0, this.sigma_x); };
  MapaMitad.prototype.fy = function (a, b) { return fraccionGaussiana(a, b, this.L / 2, this.sigma_y); };
  MapaMitad.prototype.potencia = function () {
    return this.q_pico * 2 * Math.PI * this.sigma_x * this.sigma_y *
      this.fx(0.0, this.W) * this.fy(0.0, this.L);
  };
  MapaMitad.prototype.qMedio = function () { return this.potencia() / (this.W * this.L); };
  MapaMitad.prototype.q = function (x, y) {
    return this.q_pico * Math.exp(-x * x / (2 * this.sigma_x * this.sigma_x) -
      (y - this.L / 2) * (y - this.L / 2) / (2 * this.sigma_y * this.sigma_y));
  };
  // Flujo MEDIO de cada porcion de la franja [x0, x0 + ancho], N x M, integrado
  MapaMitad.prototype.mapaNodal = function (x0, ancho, M, N) {
    var dx = ancho / M, dy = this.L / N, fx = [], fy = [], i, j;
    for (i = 0; i < M; i++) fx.push(this.fx(x0 + i * dx, x0 + (i + 1) * dx));
    for (j = 0; j < N; j++) fy.push(this.fy(j * dy, (j + 1) * dy));
    var k = this.q_pico * 2 * Math.PI * this.sigma_x * this.sigma_y / (dx * dy);
    var salida = [];
    for (j = 0; j < N; j++) {
      var fila = new Array(M);
      for (i = 0; i < M; i++) fila[i] = k * fy[j] * fx[i];
      salida.push(fila);
    }
    return salida;
  };

  // Biseccion de receptor.py
  function biseccion(f, a, b, tol, iteraciones) {
    tol = tol === undefined ? 1e-6 : tol;
    iteraciones = iteraciones === undefined ? 80 : iteraciones;
    var fa = f(a), fb = f(b);
    if (fa * fb > 0) throw new ErrorModelo('La bisección no encierra la raíz.', '');
    for (var n = 0; n < iteraciones; n++) {
      if (b - a < tol) break;
      var m = 0.5 * (a + b), fm = f(m);
      if (fa * fm <= 0) { b = m; } else { a = m; fa = fm; }
    }
    return 0.5 * (a + b);
  }

  // Pico y media dados, sigmas proporcionales a los lados del receptor
  // desarrollado: sigma_x = kappa*2W y sigma_y = kappa*L
  function mapaDesdeMedioYPico(L, w, n, q_medio, q_pico) {
    if (!(q_medio > 0 && q_medio < q_pico)) {
      throw new ErrorModelo('El flujo medio tiene que ser menor que el de pico.',
                            'Baja el flujo medio o sube el de pico.');
    }
    var W = n * w;
    var kappa = biseccion(function (k) {
      return new MapaMitad(L, w, n, q_pico, k * 2 * W, k * L).qMedio() - q_medio;
    }, 1e-3, 1e3, 1e-12);
    return new MapaMitad(L, w, n, q_pico, kappa * 2 * W, kappa * L);
  }

  // ===========================================================================
  // El receptor
  // ===========================================================================

  /* cfg:
   *   altura, D        [m]      alto de las placas y diametro del octogono
   *   circunscrito     bool     D de vertice a vertice (true) o entre caras
   *   alfa_estrella    [rad]    angulo de la estrella de 8 puntas (0 = octogono)
   *   q_medio, q_pico  [W/m2]
   *   T_ent, T_sal     [K]      entrada y salida impuestas de la sal
   *   lam              [-]      mezcla entre columnas
   *   s_t, s_2l, b_i, d_sp, delta_pp  [m]   geometria del pillow plate
   *   h_ext, eps, alfa
   */
  function resolver(cfg) {
    var nPl = N_PLACAS / 2;
    var w = cfg.circunscrito ? cfg.D * Math.sin(Math.PI / N_PLACAS)
                             : cfg.D * Math.tan(Math.PI / N_PLACAS);
    var L = cfg.altura;
    var mapa = mapaDesdeMedioYPico(L, w, nPl, cfg.q_medio, cfg.q_pico);
    var lam = cfg.lam, alfa = cfg.alfa, eps = cfg.eps, h_ext = cfg.h_ext;
    if (!(lam >= 0 && lam <= 1)) throw new ErrorModelo('λ tiene que estar entre 0 y 1.', '');
    // Estrella: cada placa se dobla en dos lados que forman alfa_e con el lado
    // del octogono. Ancho desarrollado w/cos(alfa_e) y flujo cos(alfa_e)*q''
    var alfaE = cfg.alfa_estrella || 0;
    if (!(alfaE >= 0 && alfaE <= ALFA_ESTRELLA_MAX + 1e-12)) {
      throw new ErrorModelo('El ángulo de la estrella tiene que estar entre 0 y ' +
                            (ALFA_ESTRELLA_MAX * 180 / Math.PI).toFixed(1) + '°.', '');
    }
    var cosE = Math.cos(alfaE), wpp = w / cosE;
    if (!(cfg.T_sal > cfg.T_ent + 1)) {
      throw new ErrorModelo('La salida tiene que estar por encima de la entrada.',
                            'Sube la temperatura de salida o baja la de entrada.');
    }
    verificarT(cfg.T_ent); verificarT(cfg.T_sal);

    var panel = { s_t: cfg.s_t, s_2l: cfg.s_2l, b_i: cfg.b_i, d_sp: cfg.d_sp, delta_pp: cfg.delta_pp };
    var M = Math.max(redondeoPython(wpp / panel.s_t), 1);
    var N = Math.max(redondeoPython(L / panel.s_2l), 1);
    var fracSoldadura = 2 * (Math.PI * panel.d_sp * panel.d_sp / 4) / (panel.s_t * panel.s_2l);
    var n = asignacionDeConstantes(panel.s_t, panel.s_2l, panel.d_sp, panel.b_i);
    var seccion = fChI(panel.b_i, wpp, 0.0);         // Ec. (16), w_e = 0
    var de = deI(panel.b_i);                          // Ec. (14)
    var e = panel.delta_pp;
    // Cada porcion (ancho wpp/M) se proyecta sobre un tramo w/M del lado del
    // octogono y recibe cos(alfa_e) por el flujo medio de ese tramo
    var flujos = [];
    for (var kk = 0; kk < nPl; kk++) {
      var fk = mapa.mapaNodal(kk * w, w, M, N);
      if (cosE !== 1) fk.forEach(function (fila) { for (var c = 0; c < M; c++) fila[c] *= cosE; });
      flujos.push(fk);
    }

    var dy = L / N, A = (wpp / M) * dy, A_f = A * (1.0 - fracSoldadura);

    // -- Canal interno: ec. (3) -> Re, Pr -> Nu (20) -> h ---------------------
    function hInterno(T, G) {
      var pr = propiedades(T);
      var u = G / (pr.rho * seccion);
      var Re = pr.rho * u * de / pr.mu;
      var Nu = n.n3 * Math.pow(Re, n.n4) * Math.pow(pr.Pr, n.n5);
      return { h: Nu * pr.k / de, u: u, Re: Re, Pr: pr.Pr };
    }

    // -- Pared: Newton monotono (balance decreciente y concavo en T_w) --------
    function pared(q_abs, T_f, UA) {
      var T = T_f + q_abs / UA;
      for (var it = 0; it < 50; it++) {
        var T3 = T * T * T;
        var f = q_abs - perdidas(T, A, h_ext, eps) - UA * (T - T_f);
        var df = -A * (h_ext + 4 * eps * SIGMA * T3) - UA;
        var paso = f / df;
        T -= paso;
        if (Math.abs(paso) < 1e-8) return T;
      }
      throw new ErrorModelo('La temperatura de pared no converge.', '');
    }

    // -- Una porcion: Newton con salvaguarda en la salida ----------------------
    function porcion(q_inc, T_ent, G_col, G) {
      var q_abs = alfa * q_inc, h_T_ent = entalpia(T_ent);
      var a = SAL.T_min;
      var b = Math.min(T_ent + q_abs / (G_col * cp(T_ent)), SAL.T_max);
      var T_sal = b;
      for (var it = 0; it < 30; it++) {
        var T_f = 0.5 * (T_ent + T_sal);
        var h_int = hInterno(T_f, G).h;
        var UA = A_f / (e / K_ACERO + 1.0 / h_int);
        var T_w = pared(q_abs, T_f, UA);
        var q_f = UA * (T_w - T_f);
        var dP = A * (h_ext + 4 * eps * SIGMA * T_w * T_w * T_w);
        var D = q_f - G_col * (entalpia(T_sal) - h_T_ent);
        if (it === 0 && D > 0) {
          throw new ErrorModelo(
            'La sal se saldría de ' + (SAL.T_max - CERO_CELSIUS).toFixed(0) + ' °C en una porción.',
            'Baja la temperatura de salida o el flujo incidente.');
        }
        if (D > 0) a = T_sal; else b = T_sal;
        var dD = -G_col * cp(T_sal) - 0.5 * UA * dP / (dP + UA);
        var paso = -D / dD;
        if (Math.abs(paso) < TOL_FILA) {
          return { T_sal: T_sal, T_w: T_w, T_pel: T_f + q_f / (A_f * h_int), q: q_f };
        }
        var nuevo = T_sal + paso;
        T_sal = (nuevo <= a || nuevo >= b) ? 0.5 * (a + b) : nuevo;
      }
      throw new ErrorModelo('La temperatura de salida de una porción no converge.', '');
    }

    // -- Mezcla de la fila con lambda, en entalpia (columnas de igual gasto) --
    function mezclar(T_ad) {
      if (lam === 1.0) return T_ad.slice();
      var h = new Array(M), media = 0.0, i;
      for (i = 0; i < M; i++) { h[i] = entalpia(T_ad[i]); media += h[i]; }
      media /= M;
      var salida = new Array(M);
      for (i = 0; i < M; i++) salida[i] = TdesdeH(lam * h[i] + (1.0 - lam) * media);
      return salida;
    }

    // -- Una placa ---------------------------------------------------------------
    function placa(k, T_ent, G) {
      var G_col = G / M, sube = (k % 2 === 0);
      var flujo = flujos[k];
      var T_fluido = new Array(N), T_pared = new Array(N), T_pel = new Array(N);
      var T_adiab = new Array(N), entradas = new Array(N), q_nodo = new Array(N);
      var Q_util = 0, Q_conv = 0, Q_rad = 0, i, s, j;
      var Tc4 = Math.pow(T_CIELO, 4);
      var T = new Array(M);
      for (i = 0; i < M; i++) T[i] = T_ent;
      for (s = 0; s < N; s++) {
        j = sube ? s : N - 1 - s;
        entradas[j] = T.slice();
        var T_ad = new Array(M), Tw = new Array(M), Tp = new Array(M), qf = new Array(M);
        for (i = 0; i < M; i++) {
          var r = porcion(flujo[j][i] * A, T[i], G_col, G);
          T_ad[i] = r.T_sal; Tw[i] = r.T_w; Tp[i] = r.T_pel; qf[i] = r.q;
          Q_util += r.q;
          Q_conv += h_ext * (r.T_w - T_AMB) * A;
          Q_rad += eps * SIGMA * (Math.pow(r.T_w, 4) - Tc4) * A;
        }
        T = mezclar(T_ad);
        T_fluido[j] = T; T_pared[j] = Tw; T_pel[j] = Tp; T_adiab[j] = T_ad; q_nodo[j] = qf;
      }
      // Colector: mezcla completa, sea cual sea lambda
      var hm = 0;
      for (i = 0; i < M; i++) hm += entalpia(T[i]);
      var T_sal = TdesdeH(hm / M);
      var Q_inc = 0;
      for (j = 0; j < N; j++) for (i = 0; i < M; i++) Q_inc += flujo[j][i] * A;

      // -- Perdida de carga: ec. (19) fila a fila, mas la hidrostatica ---------
      var G2 = Math.pow(G / seccion, 2);
      var fr = 0, di = 0, ac = 0, hi = 0;
      // Validez: Re y Pr de cada porcion frente a los rangos de las correlaciones
      var val = { Re: [Infinity, -Infinity], Pr: [Infinity, -Infinity],
                  fueraPiperRe: 0, fueraPiperPr: 0, fueraArsenyevaRe: 0 };
      for (i = 0; i < M; i++) {
        var frI = 0, hiI = 0;
        for (j = 0; j < N; j++) {
          var pr = propiedades(0.5 * (entradas[j][i] + T_fluido[j][i]));
          var u = G / (pr.rho * seccion);
          var Re = pr.rho * u * de / pr.mu;
          frI += n.n1 * Math.pow(Re, n.n2) * dy / de * pr.rho * u * u / 2;
          hiI += pr.rho * dy;
          if (Re < val.Re[0]) val.Re[0] = Re;
          if (Re > val.Re[1]) val.Re[1] = Re;
          if (pr.Pr < val.Pr[0]) val.Pr[0] = pr.Pr;
          if (pr.Pr > val.Pr[1]) val.Pr[1] = pr.Pr;
          if (Re < RANGOS.piper.Re[0] || Re > RANGOS.piper.Re[1]) val.fueraPiperRe++;
          if (pr.Pr < RANGOS.piper.Pr[0] || pr.Pr > RANGOS.piper.Pr[1]) val.fueraPiperPr++;
          if (Re < RANGOS.arsenyeva.Re[0] || Re > RANGOS.arsenyeva.Re[1]) val.fueraArsenyevaRe++;
        }
        var v_ent = 1.0 / rho(entradas[sube ? 0 : N - 1][i]);
        var v_sal = 1.0 / rho(T[i]);
        fr += frI;
        di += ZETA_DZ / 2 * G2 * (v_ent + v_sal);
        ac += G2 * (v_sal - v_ent);
        hi += (sube ? 1 : -1) * GRAVEDAD * hiI;
      }
      var perfil = T.slice();
      return {
        numero: k + 1, sube: sube, x0: k * wpp, T_ent: T_ent, T_sal: T_sal,
        perfilSalida: perfil, T_fluido: T_fluido, T_pared: T_pared, T_pelicula: T_pel,
        T_adiabatica: T_adiab, entradas: entradas, q_nodo: q_nodo, flujo: flujo,
        Q_inc: Q_inc, Q_util: Q_util, Q_conv: Q_conv, Q_rad: Q_rad, Q_refl: (1 - alfa) * Q_inc,
        dp_friccion: fr / M, dp_distribucion: di / M, dp_aceleracion: ac / M,
        dp_hidrostatica: hi / M, validez: val
      };
    }

    function recorrer(G) {
      var placas = [], T = cfg.T_ent;
      for (var k = 0; k < nPl; k++) { placas.push(placa(k, T, G)); T = placas[k].T_sal; }
      return placas;
    }

    // -- Gasto inicial con las perdidas estimadas (gasto_inicial) --------------
    var dh = entalpia(cfg.T_sal) - entalpia(cfg.T_ent);
    function gastoInicial() {
      var h_ent = entalpia(cfg.T_ent);
      var G = alfa * mapa.potencia() / dh;
      var Q = [], sumaQ = 0, k, i, j;
      for (k = 0; k < nPl; k++) {
        var s = 0;
        for (j = 0; j < N; j++) for (i = 0; i < M; i++) s += flujos[k][j][i] * A;
        Q.push(s); sumaQ += s;
      }
      var perd = 0, acumulado = 0;
      for (k = 0; k < nPl; k++) {
        var T_f = TdesdeH(h_ent + dh * (acumulado + Q[k] / 2) / sumaQ);
        acumulado += Q[k];
        var R = e / K_ACERO + 1.0 / hInterno(T_f, G).h;
        for (j = 0; j < N; j++) for (i = 0; i < M; i++) {
          var T_w = T_f + alfa * flujos[k][j][i] * R / (1.0 - fracSoldadura);
          perd += perdidas(T_w, A, h_ext, eps);
        }
      }
      return (alfa * sumaQ - perd) / dh;
    }

    // -- Punto fijo del gasto, acelerado con la secante --------------------------
    var G = gastoInicial(), anterior = null, placas = null, iteraciones = 0;
    for (var it = 1; it <= 30; it++) {
      placas = recorrer(G);
      if (Math.abs(placas[nPl - 1].T_sal - cfg.T_sal) < 1e-3) { iteraciones = it; break; }
      var Qu = 0;
      for (var q = 0; q < nPl; q++) Qu += placas[q].Q_util;
      var r = Qu / dh - G, Gn;
      if (anterior === null || r === anterior[1]) Gn = G + r;
      else Gn = G - r * (G - anterior[0]) / (r - anterior[1]);
      anterior = [G, r]; G = Gn;
      if (!(G > 0)) throw new ErrorModelo('El gasto no converge.', 'Revisa el punto de trabajo.');
    }
    if (!iteraciones) throw new ErrorModelo('El gasto no converge en 30 pasadas.', '');

    // -- Resumen ----------------------------------------------------------------
    var res = {
      mapa: mapa, w: w, wpp: wpp, alfaEstrella: alfaE, L: L, W: nPl * wpp, nPlacas: nPl,
      M: M, N: N, dy: dy, dx: wpp / M,
      A_celda: A, A_f: A_f, fracSoldadura: fracSoldadura, n: n, panel: panel,
      seccion: seccion, de: de, lam: lam, alfa: alfa, G: G, iteraciones: iteraciones,
      placas: placas, T_ent: cfg.T_ent, T_sal: placas[nPl - 1].T_sal
    };
    var s = { Q_inc: 0, Q_util: 0, Q_conv: 0, Q_rad: 0, Q_refl: 0, dp_friccion: 0,
              dp_distribucion: 0, dp_aceleracion: 0, dp_hidrostatica: 0 };
    placas.forEach(function (p) { for (var c in s) s[c] += p[c]; });
    for (var c in s) res[c] = s[c];
    res.dp = s.dp_friccion + s.dp_distribucion + s.dp_aceleracion;
    res.rendimiento = s.Q_util / s.Q_inc;
    res.entrada = hInterno(cfg.T_ent, G);
    res.salida = hInterno(cfg.T_sal, G);
    res.rhoMedia = 0.5 * (rho(cfg.T_ent) + rho(cfg.T_sal));
    res.potenciaBombeo = G * res.dp / res.rhoMedia;              // de la mitad [W]
    res.validez = validez(placas, M * N * nPl);

    // Campos del receptor desarrollado, N filas x (nPl*M) columnas
    var campos = { flujo: [], sal: [], pared: [], pelicula: [] }, j, i, k2;
    var maxPared = -Infinity, maxPel = -Infinity;
    for (j = 0; j < N; j++) {
      var fq = [], fs = [], fp = [], fl = [];
      for (k2 = 0; k2 < nPl; k2++) {
        var p = placas[k2];
        for (i = 0; i < M; i++) {
          fq.push(p.flujo[j][i]); fs.push(p.T_fluido[j][i]);
          fp.push(p.T_pared[j][i]); fl.push(p.T_pelicula[j][i]);
          if (p.T_pared[j][i] > maxPared) maxPared = p.T_pared[j][i];
          if (p.T_pelicula[j][i] > maxPel) maxPel = p.T_pelicula[j][i];
        }
      }
      campos.flujo.push(fq); campos.sal.push(fs); campos.pared.push(fp); campos.pelicula.push(fl);
    }
    res.campos = campos;
    res.T_pared_max = maxPared;
    res.T_pelicula_max = maxPel;
    res.limiteAcero = maxPared > LIMITE_AISI321;
    res.limiteSal = maxPel > T_PELICULA_MAX;
    res.xCentros = []; res.yCentros = [];
    for (i = 0; i < nPl * M; i++) res.xCentros.push((i + 0.5) * wpp / M);
    for (j = 0; j < N; j++) res.yCentros.push((j + 0.5) * dy);

    // Recorrido de la sal, fila a fila en el sentido del fluido
    var recorrido = { s: [], sal: [], pared: [], pelicula: [] };
    placas.forEach(function (p, kp) {
      for (var st = 0; st < N; st++) {
        var jj = p.sube ? st : N - 1 - st, ms = 0, mp = -Infinity, ml = -Infinity;
        for (var ii = 0; ii < M; ii++) {
          ms += p.T_fluido[jj][ii];
          if (p.T_pared[jj][ii] > mp) mp = p.T_pared[jj][ii];
          if (p.T_pelicula[jj][ii] > ml) ml = p.T_pelicula[jj][ii];
        }
        recorrido.s.push(kp * L + (st + 0.5) * dy);
        recorrido.sal.push(ms / M); recorrido.pared.push(mp); recorrido.pelicula.push(ml);
      }
    });
    res.recorrido = recorrido;
    return res;
  }

  // ===========================================================================
  // Validez de las correlaciones
  // ===========================================================================

  // Junta lo de las cuatro placas y redacta un aviso por cada rango que se
  // incumple. No para el calculo: el modelo sigue dando numeros, pero con Nu y
  // f extrapolados, y eso hay que saberlo al leerlos.
  function validez(placas, total) {
    var v = { Re: [Infinity, -Infinity], Pr: [Infinity, -Infinity],
              fueraPiperRe: 0, fueraPiperPr: 0, fueraArsenyevaRe: 0 };
    placas.forEach(function (p) {
      var w = p.validez;
      v.Re[0] = Math.min(v.Re[0], w.Re[0]); v.Re[1] = Math.max(v.Re[1], w.Re[1]);
      v.Pr[0] = Math.min(v.Pr[0], w.Pr[0]); v.Pr[1] = Math.max(v.Pr[1], w.Pr[1]);
      v.fueraPiperRe += w.fueraPiperRe; v.fueraPiperPr += w.fueraPiperPr;
      v.fueraArsenyevaRe += w.fueraArsenyevaRe;
    });
    v.fueraPiperRe /= total; v.fueraPiperPr /= total; v.fueraArsenyevaRe /= total;

    var avisos = [];
    function sentido(valor, rango) {
      return valor[0] > rango[1] ? 'por encima' : valor[1] < rango[0] ? 'por debajo' : 'en parte fuera';
    }
    if (v.fueraPiperRe > 0) {
      avisos.push({ clave: 'piperRe', grave: true, fuente: RANGOS.piper.nombre,
        magnitud: 'Re', rango: RANGOS.piper.Re, valor: v.Re, fraccion: v.fueraPiperRe,
        sentido: sentido(v.Re, RANGOS.piper.Re), ecuaciones: RANGOS.piper.ecuaciones });
    }
    if (v.fueraPiperPr > 0) {
      avisos.push({ clave: 'piperPr', grave: true, fuente: RANGOS.piper.nombre,
        magnitud: 'Pr', rango: RANGOS.piper.Pr, valor: v.Pr, fraccion: v.fueraPiperPr,
        sentido: sentido(v.Pr, RANGOS.piper.Pr), ecuaciones: RANGOS.piper.ecuaciones });
    }
    if (v.fueraArsenyevaRe > 0) {
      avisos.push({ clave: 'arsenyevaRe', grave: false, fuente: RANGOS.arsenyeva.nombre,
        magnitud: 'Re', rango: RANGOS.arsenyeva.Re, valor: v.Re, fraccion: v.fueraArsenyevaRe,
        sentido: sentido(v.Re, RANGOS.arsenyeva.Re), ecuaciones: RANGOS.arsenyeva.ecuaciones });
    }
    v.avisos = avisos;
    return v;
  }

  var ROC = {
    RANGOS: RANGOS,
    CERO_CELSIUS: CERO_CELSIUS, SIGMA: SIGMA, T_AMB: T_AMB, T_CIELO: T_CIELO,
    LIMITE_AISI321: LIMITE_AISI321, T_PELICULA_MAX: T_PELICULA_MAX, K_ACERO: K_ACERO,
    ZETA_DZ: ZETA_DZ, N_PLACAS: N_PLACAS, ALFA_ESTRELLA_MAX: ALFA_ESTRELLA_MAX, LAMBDA: 0.7, SAL: SAL,
    ErrorModelo: ErrorModelo, FAMILIAS: FAMILIAS,
    asignacionDeConstantes: asignacionDeConstantes, propiedades: propiedades,
    entalpia: entalpia, TdesdeH: TdesdeH, mapaDesdeMedioYPico: mapaDesdeMedioYPico,
    resolver: resolver
  };

  global.ROC = ROC;
  if (typeof module !== 'undefined' && module.exports) module.exports = ROC;
})(typeof globalThis !== 'undefined' ? globalThis : this);
