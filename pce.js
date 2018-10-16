
var generators = {
  range: (n) => [...Array(n).keys()],
  cartesian: (a, b, ...c) => (b ? generators.cartesian([].concat(...a.map(d => b.map(e => [].concat(d, e)))), ...c) : a)
};

var reductions = {
  sum: (a, b) => a + b,
  product: (a, b) => a * b
};

const factorial = (n) => (n == 0 || n == 1) ? 1 : factorial(n-1) * n;

const germType = { NORMAL: 0, UNIFORM: 1 };

class PCExpansion {

  constructor(f, germs, orders, totalOrder) {
    this.f = f;
    this.germs = germs;
    this.orders = orders;
    this.totalOrder = totalOrder;
    this.multiindexSet = (totalOrder > 0) 
      ? generators.cartesian.apply(null, this.orders.map(order => generators.range(order+1))).filter(multiindex => multiindex.reduce(reductions.sum) <= totalOrder)
      : generators.cartesian.apply(null, this.orders.map(order => generators.range(order+1)));
    this.quads = this.germs.map((germ, i) => generators.range(this.orders[i] + 3).map(j => new GaussianQuadrature(germ, j + 1)));
    this.tensorQuads = this.multiindexSet.map((multiindex) => ({ 
      abscissae: generators.cartesian.apply(null, this.germs.map((germ, i) => this.quads[i][multiindex[i]+2].abscissae)), 
      weights: generators.cartesian.apply(null, this.germs.map((germ, i) => this.quads[i][multiindex[i]+2].weights)) 
    }));
    this.norms = this.multiindexSet.map((multiindex, k) => this.germs.map((germ, i) => OrthogonalPolynomial.norm(germ, multiindex[i])).reduce(reductions.product));
    this.coeffs = this.multiindexSet.map((multiindex, k) => 
      this.tensorQuads[k].abscissae.map((abscissa, j) => 
        this.tensorQuads[k].weights[j].reduce(reductions.product) * this.f(abscissa) * this.evalTensorPoly(multiindex, abscissa)
      ).reduce(reductions.sum) / this.norms[k]
    );
  }

  evalTensorPoly(multiindex, x) {
    return this.germs.map((germ, i) => OrthogonalPolynomial.evaluate(germ, multiindex[i], x[i])).reduce(reductions.product);
  }

  evaluate(x) {
    return this.coeffs.map((coeff, k) => coeff * this.evalTensorPoly(this.multiindexSet[k], x)).reduce(reductions.sum);
  }

  getMoments() {
    return { 
      mean: this.coeffs[0], 
      variance: generators.range(this.coeffs.length).slice(1).map(k => Math.pow(this.coeffs[k], 2) * this.norms[k]).reduce(reductions.sum) 
    };
  }

  getTensorPolyLatex(multiindex) {
    return this.germs.map(function (type, i) {
      switch (type) {
        case germType.NORMAL: return 'H_{' + multiindex[i] + '}(\\xi_{' + (i+1) + '})';
        case germType.UNIFORM: return 'P_{' + multiindex[i] + '}(\\xi_{' + (i+1) + '})';
      }
    }).join('');
  }

  getGermsLatex() {
    return this.germs.map(function (type, i) {
      switch (type) {
        case germType.NORMAL: return '\\xi_{' + (i+1) + '}\\sim\\mathcal{N}(0,1)';
        case germType.UNIFORM: return '\\xi_{' + (i+1) + '}\\sim\\mathcal{U}(-1,1)';
      }
    }).join(';\\\\');
  }

  getOrdersLatex() {
    return this.orders.map((order, i) => '0\\leq j_{' + (i+1) + '}\\leq ' + order).join('\\text{ and }'); 
  }

}

class OrthogonalPolynomial {

  static evaluate(type, order, x) {
    if (order == 0) return 1;
    if (order == 1) return x;
    let yn2 = 1, yn1 = x, y = 0;
    for (let k = 2; k <= order; k++) {
      switch (type) {
        case germType.NORMAL: y = x * yn1 - (k - 1) * yn2; break; // Hermite (probabilist's)
        case germType.UNIFORM: y = ((2 * k - 1) * x * yn1 - (k - 1) * yn2) / k; break; // Legendre
      }
      yn2 = yn1;
      yn1 = y;
    }
    return y;
  }

  static norm(type, order) {
    switch (type) {
      case germType.NORMAL: return factorial(order); break; // Hermite
      case germType.UNIFORM: return 2 / (2 * order + 1); break; // Legendre
    }
  }

}

class GaussianQuadrature {

  constructor(type, order) {
    this.order = order;
    switch (type) {
      case germType.NORMAL: 
        this.abscissae = GaussianQuadrature.hermite[order - 1].abscissae;
        this.weights = GaussianQuadrature.hermite[order - 1].weights;
        break;
      case germType.UNIFORM:
        this.abscissae = GaussianQuadrature.legendre[order - 1].abscissae;
        this.weights = GaussianQuadrature.legendre[order - 1].weights;
        break;
    }
  }

}

GaussianQuadrature.hermite = [
  { abscissae: [0], weights: [1]},
  { abscissae: [-1, 1], weights: [0.5, 0.5]},
  { abscissae: [-1.73205080756888, -2.06976446027777e-16, 1.73205080756888], weights: [0.166666666666667, 0.666666666666667, 0.166666666666667]},
  { abscissae: [-2.33441421833898, -0.741963784302726, 0.741963784302726, 2.33441421833898], weights: [0.0458758547680685, 0.454124145231932, 0.454124145231932, 0.0458758547680686]},
  { abscissae: [-2.85697001387281, -1.35562617997427, 3.86574996532659e-17, 1.35562617997427, 2.85697001387281], weights: [0.0112574113277207, 0.222075922005613, 0.533333333333333, 0.222075922005613, 0.0112574113277207]},
  { abscissae: [-3.32425743355212, -1.88917587775371, -0.616706590192594, 0.616706590192594, 1.88917587775371, 3.32425743355212], weights: [0.00255578440205626, 0.0886157460419144, 0.408828469556029, 0.408828469556029, 0.0886157460419144, 0.00255578440205624]},
  { abscissae: [-3.75043971772574, -2.36675941073454, -1.15440539473997, 1.15327337087898e-16, 1.15440539473997, 2.36675941073454, 3.75043971772574], weights: [0.000548268855972218, 0.0307571239675865, 0.240123178605013, 0.457142857142857, 0.240123178605013, 0.0307571239675865, 0.000548268855972215]},
  { abscissae: [-4.14454718612589, -2.80248586128754, -1.63651904243511, -0.539079811351375, 0.539079811351375, 1.63651904243511, 2.80248586128754, 4.14454718612589], weights: [0.000112614538375368, 0.00963522012078822, 0.117239907661759, 0.373012257679077, 0.373012257679077, 0.117239907661759, 0.00963522012078826, 0.000112614538375368]},
  { abscissae: [-4.51274586339979, -3.20542900285647, -2.07684797867783, -1.02325566378913, -5.69133976269097e-17, 1.02325566378913, 2.07684797867783, 3.20542900285647, 4.51274586339979], weights: [2.23458440077466e-05, 0.00278914132123177, 0.0499164067652179, 0.244097502894939, 0.406349206349206, 0.244097502894939, 0.049916406765218, 0.00278914132123176, 2.23458440077465e-05]},
  { abscissae: [-4.85946282833231, -3.58182348355193, -2.48432584163895, -1.46598909439116, -0.484935707515498, 0.484935707515498, 1.46598909439116, 2.48432584163896, 3.58182348355193, 4.85946282833231], weights: [4.31065263071827e-06, 0.000758070934312221, 0.0191115805007703, 0.135483702980268, 0.344642334932019, 0.344642334932019, 0.135483702980268, 0.0191115805007703, 0.00075807093431222, 4.3106526307183e-06]}
];

GaussianQuadrature.legendre = [
  { abscissae: [0], weights: [2]},
  { abscissae: [-0.577350269189626, 0.577350269189626], weights: [1, 1]},
  { abscissae: [-0.774596669241484, 0, 0.774596669241484], weights: [0.555555555555555, 0.888888888888889, 0.555555555555555]},
  { abscissae: [-0.861136311594053, -0.339981043584856, 0.339981043584856, 0.861136311594053], weights: [0.347854845137454, 0.652145154862546, 0.652145154862546, 0.347854845137454]},
  { abscissae: [-0.906179845938664, -0.538469310105683, 0, 0.538469310105683, 0.906179845938664], weights: [0.236926885056189, 0.478628670499366, 0.568888888888889, 0.478628670499366, 0.236926885056189]},
  { abscissae: [-0.932469514203152, -0.661209386466264, -0.238619186083197, 0.238619186083197, 0.661209386466264, 0.932469514203152], weights: [0.17132449237917, 0.360761573048139, 0.467913934572691, 0.467913934572691, 0.360761573048139, 0.17132449237917]},
  { abscissae: [-0.949107912342759, -0.741531185599394, -0.405845151377397, 0, 0.405845151377397, 0.741531185599394, 0.949107912342759], weights: [0.129484966168869, 0.279705391489277, 0.381830050505119, 0.417959183673469, 0.381830050505119, 0.279705391489277, 0.129484966168869]},
  { abscissae: [-0.960289856497536, -0.796666477413627, -0.525532409916329, -0.18343464249565, 0.18343464249565, 0.525532409916329, 0.796666477413627, 0.960289856497536], weights: [0.101228536290376, 0.222381034453374, 0.313706645877887, 0.362683783378362, 0.362683783378362, 0.313706645877887, 0.222381034453374, 0.101228536290376]},
  { abscissae: [-0.968160239507626, -0.836031107326636, -0.61337143270059, -0.324253423403809, 0, 0.324253423403809, 0.61337143270059, 0.836031107326636, 0.968160239507626], weights: [0.0812743883615746, 0.180648160694858, 0.260610696402935, 0.312347077040003, 0.33023935500126, 0.312347077040003, 0.260610696402935, 0.180648160694858, 0.0812743883615746]},
  { abscissae: [-0.973906528517172, -0.865063366688985, -0.679409568299024, -0.433395394129247, -0.148874338981631, 0.148874338981631, 0.433395394129247, 0.679409568299024, 0.865063366688985, 0.973906528517172], weights: [0.0666713443086879, 0.149451349150581, 0.219086362515982, 0.269266719309996, 0.295524224714753, 0.295524224714753, 0.269266719309996, 0.219086362515982, 0.149451349150581, 0.0666713443086879]}
];
