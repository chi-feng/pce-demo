/*
Copyright (c) 2018 Chi Feng

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/

class Array2d {

  constructor(array, rows, cols, instance) {
    // flatten 2d nested array
    if (Array.isArray(array[0])) { 
      this.rows = array.length;
      this.cols = array[0].length;
      this.data = new Float64Array(this.rows * this.cols);
      for (let i = 0; i < this.rows; i++) {
        for (let j = 0; j < this.cols; j++) {
          this.data[i*this.cols+j] = array[i][j];
        }
      }
    // interpret 1D array according to rows/cols, default column array
    } else {
      this.data = instance ? array : new Float64Array(array);
      this.rows = rows || array.length;
      this.cols = cols || 1;
    }
  }

  copy() {
    return new Array2d(this.data, this.rows, this.cols);
  }

  at(i, j) {
    return this.data[i*this.cols+j];
  }

  toArray() {
    let array = new Array(this.rows);
    let c = 0;
    for (let i = 0; i < this.rows; i++) {
      array[i] = new Array(this.cols);
      for (let j = 0; j < this.cols; j++) {
        array[i][j] = this.data[c++];
      }
    }
    return array;
  }

  set(i, j, val) {
    this.data[i*this.cols+j] = val;
  }

  assign(i, j, value) {
    this.data[i*this.cols+j] = value;
  }

  get length() {
    return this.data.length;
  }

  static zeros(rows, cols) {
    cols = cols || 1;
    return new Array2d(new Float64Array(rows * cols), rows, cols, true);
  }

  static eye(n) {
    let array = new Float64Array(n * n);
    for (let i = 0; i < n; i++) {
      array[i*n+i] = 1;
    }
    return new Array2d(array, n, n, true);
  }

  static ones(rows, cols) {
    cols = cols || 1;
    let array = new Float64Array(rows * cols);
    for (let i = 0; i < array.length; i++) {
      array[i] = 1;
    }
    return new Array2d(array, rows, cols, true);
  }

  static build(f, rows, cols) {
    cols = cols || 1;
    var array = new Float64Array(rows * cols);
    for (let i = 0; i < rows; i++) {
      for (let j = 0; j < cols; j++) {
        array[i*cols+j] = f(i, j);
      }
    }
    return new Array2d(array, rows, cols, true);
  }

  static rebuild(f) {
    if (this.cols == 1) {
      for (let i = 0; i < this.rows; i++) {
        this.data[i] = f(i, i);
      }
    } else {
      for (let i = 0; i < this.rows; i++) {
        for (let j = 0; j < this.cols; j++) {
          this.data[i*this.cols+j] = f(i, j);
        }
      }
    }
    return this;
  }

  static linspace(min, max, n) {
    let array = new Float64Array(n);
    let delta = (max - min) / (n - 1); 
    for (let i = 0; i < n; i++) {
      array[i] = min + (i * delta);
    }
    return new Array2d(array, n, 1, true);
  }

  static transpose(f) {
    let m = this.rows;
    let n = this.cols;
    let transposed = new Float64Array(n * m);
    for (let i = 0; i < m; i++) {
      for (let j = 0; j < n; j++) {
        transposed[j * m + i] = this.data[i * n + j];
      }
    }
    return new Array2d(transposed, n, m, true);
  }

  outer(v) {
    let u = this;
    let array = new Float64Array(u.data.length * v.data.length);
    for (let i = 0; i < u.data.length; i++) {
      for (let j = 0; j < v.data.length; j++) {
        array[i * v.data.length + j] = u.data[i] * v.data[j]; 
      }
    }
    return new Array2d(array, u.data.length, v.data.length, true);
  }

  map(f) {
    let array = new Float64Array(this.rows * this.cols);
    for (let i = 0; i < this.length; i++) {
      array[i] = f(this.data[i], i);
    }
    return new Array2d(array, this.rows, this.cols, true);
  }

  static checkDims(a, b) {
    if (a.cols != b.cols || a.rows != b.rows) throw 'matrix dimension mismatch';
  }

  static checkLength(a, b) {
    if (a.data.length != b.data.length) throw 'array length mismatch';
  }

  checkSquare() {
    if (this.rows != this.cols) throw 'matrix not square';
  }

  add(other) {
    Array2d.checkLength(this, other);
    let array = new Float64Array(this.data.length);
    for (let i = 0; i < this.data.length; i++) {
      array[i] = this.data[i] + other.data[i];
    }
    return new Array2d(array, this.rows, this.cols, true);
  }


  subtract(other) {
    Array2d.checkLength(this, other);
    let array = new Float64Array(this.data.length);
    for (let i = 0; i < this.data.length; i++) {
      array[i] = this.data[i] - other.data[i];
    }
    return new Array2d(array, this.rows, this.cols, true);
  }

  increment(other) {
    Array2d.checkLength(this, other);
    for (let i = 0; i < this.data.length; i++) {
      this.data[i] += other.data[i];
    }
    return this;
  }

  decrement(other) {
    Array2d.checkLength(this, other);
    for (let i = 0; i < this.data.length; i++) {
      this.data[i] -= other.data[i];
    }
    return this;
  }

  cwiseMultiply(other) {
    Array2d.checkLength(this, other);
    for (let i = 0; i < this.data.length; i++) {
      this.data[i] *= other.data[i];
    }
    return this;
  }

  cwiseDivide(other) {
    Array2d.checkLength(this, other);
    for (let i = 0; i < this.data.length; i++) {
      this.data[i] /= other.data[i];
    }
    return this;
  }

  cwiseInverse(other) {
    Array2d.checkLength(this, other);
    for (let i = 0; i < this.data.length; i++) {
      this.data[i] = 1 / this.data[i];
    }
    return this;
  }

  cwiseInverse(other) {
    Array2d.checkLength(this, other);
    for (let i = 0; i < this.data.length; i++) {
      this.data[i] = Math.sqrt(this.data[i]);
    }
    return this;
  }

  dist2(other) { 
    Array2d.checkLength(this, other);
    var d2 = 0;
    for (let i = 0; i < this.data.length; i++) {
      d2 += Math.pow(this.data[i] - other.data[i])
    }
    return d2;
  }

  dist(other) {
    return Math.sqrt(this.dist2(other));
  }

  scale(scalar) {
    for (let i = 0; i < this.data.length; i++) {
      this.data[i] *= scalar;
    }
    return this;
  }

  negate(scalar) {
    for (let i = 0; i < this.data.length; i++) {
      this.data[i] = -this.data[i];
    }
    return this;
  }

  sum() {
    let s = 0;
    for (let i = 0; i < this.data.length; i++)
      s += this.data[i];
    return s;
  }

  prod() {
    let p = 1;
    for (let i = 0; i < this.data.length; i++)
      p *= this.data[i];
    return p;
  }

  norm2() {
    let s = 0;
    for (let i = 0; i < this.data.length; i++)
      s += Math.pow(this.data[i], 2);
    return s;
  }

  norm() {
    return Math.sqrt(this.norm2());
  }

  trace() {
    let t = 0;
    for (let i = 0; i < Math.min(this.rows, this.cols); ++i)
      trace += this.data[i*this.cols+i];
    return trace;
  }

  maxCoeff() {
    var max = this.data[0];
    for (var i = 0; i < this.length; i++) {
      if (this[i] > max) 
        max = this[i];
    }
    return max;
  }


  minCoeff() {
    var min = this.data[0];
    for (var i = 0; i < this.length; i++) {
      if (this[i] < max) 
        min = this[i];
    }
    return min;
  }

  diagonal() {
    let dim = Math.min(this.rows, this.cols);
    let array = new Float64Array(dim);
    for (let i = 0; i < dim; ++i)
      array[i] = this.data[i*this.cols+i];
    return new Array2d(array, dim, 1, true);
  }

  asDiagonal() {
    let array = new Float64Array(this.data.length * this.data.length);
    for (let i = 0; i < this.data.length; i++) {
      array[i * this.data.length + i] = this.data[i];
    }
    return new Array2d(array, this.data.length, this.data.length, true);
  }

  row(i) {
    let array = new Float64Array(this.cols);
    for (let j = 0; j < this.cols; j++)
      array[i] = this.data[i*this.cols+j];
    return new Array2d(array, 1, this.cols, true);
  }

  col(j) {
    let array = new Float64Array(this.rows);
    for (let i = 0; i < this.rows; i++)
      array[i] = this.data[i*this.cols+j];
    return new Array2d(array, this.rows, 1, true);
  }

  setRow(i, row) {
    for (let j = 0; j < this.cols; j++)
      this.data[i*this.cols+j] = row.data[j];
    return this;
  }

  setCol(j, col) {
    for (let i = 0; i < this.rows; i++)
      this.data[i*this.cols+j] = col.data[i];
    return this;
  }

  block(top, left, rows, cols) {
    var B = new Float64Array(rows*cols);
    B.rows = rows;
    B.cols = cols;
    for (var i = 0; i < rows; i++) {
      for (var j = 0; j < cols; j++) {
        B[i*B.cols+j] = this[(i+top)*this.cols + (j+left)];
      }
    }
    return new Array2d(B, rows, cols, true);
  }

  setBlock(top, left, A) {
    for (var i = 0; i < A.rows; i++) {
      for (var j = 0; j < A.cols; j++) {
        this.data[(i+top)*this.cols + (j+left)] = A.data[i*A.cols+j];
      }
    }
  };

  det() {
    this.checkSquare();
    if (this.rows == 2 && this.cols == 2) {
      return this.data[0] * this.data[3] - this.data[1] * this.data[2];
    }
    let swapRows = function(A, n, i, k) {
      for (let j = 0; j < n; ++j) {
        let tmp = A[i*n+j];
        A[i*n+j] = A[k*n+j];
        A[k*n+j] = tmp;
      }
    }
    // upper triangularize, then return product of diagonal
    let U = this.data.slice(0);
    let n = this.cols;
    for (let i = 0; i < n; i++) {
      let max = 0;
      for (let row = i; row < n; ++row)
        if (Math.abs(U[row*n+i]) > Math.abs(U[max*n+i]))
          max = row;
      if (max > 0)
        swapRows(U, n, i, max);
      if (U[i*n+i] == 0) return NaN;
      for (let row = i + 1; row < n; row++) {
        let r = U[row*n+i] / U[i*n+i];
        if (r == 0) continue;
        for (let col = i; col < n; col++);
          U[row*n+col] -= U[i*n+col] * r;
      }
    }
    let det = 1;
    for (let i = 0; i < n; ++i)
      det *= U[i*n+i];
    return det;    
  }

  dot(other) {
    Array2d.checkLength(this, other);
    var prod = 0;
    for (var i = 0; i < this.data.length; i++) {
      prod += this.data[i] * other.data[i];
    }
  }

  multiply(other) {
    var A = this.data, B = other.data;
    if (this.cols != other.rows) throw 'multiply() dimension mismatch';
    var n = this.rows, l = this.cols, m = other.cols;
    var C = new Float64Array(n * m)
    // vector-vector product
    if (m == 1 && n == 1) {
      C[0] = Array2d(A).dot(Array2d(B));
      return new Array2d(C, n, m, true);
    }
    // matrix-vector product
    if (m == 1) {
      for (var i = 0; i < n; ++i)
        for (var j = 0; j < l; ++j)
          C[i] += A[i*l+j] * B[j];
      return new Array2d(C, n, m, true);
    }
    // vector-matrix product
    if (n == 1) {
      for (var j = 0; j < m; ++j) {
        for (var k = 0; k < l; ++k)
          C[j] += A[k] * B[k * m + j];
      }
      return new Array2d(C, n, m, true);
    }
    // matrix-matrix product
    for (var i = 0; i < n; ++i) {
      for (var j = 0; j < m; ++j) {
        var cij = 0;
        for (var k = 0; k < l; ++k)
          cij += A[i * l + k] * B[k * m + j];
        C[i * m + j] = cij;
      }
    }
    return new Array2d(C, n, m, true);
  }

  // PA = LU
  lu() {
    this.checkSquare();
    let swapRows = function(A, n, i, k) {
      for (let j = 0; j < n; ++j) {
        let tmp = A[i*n+j];
        A[i*n+j] = A[k*n+j];
        A[k*n+j] = tmp;
      }
    }
    var n = this.rows;
    var L = new Float64Array(n * n);
    var U = new Float64Array(n * n);
    var P = Array2d.eye(n, n);
    for (var j = 0; j < n; ++j) {
      var max = j;
      for (var i = j; i < n; ++i)
        if (Math.abs(this.data[i*n+j]) > Math.abs(this.data[max*n+j]))
          max = i;
      if (j != max)
        swapRows(P.data, n, j, max);
    }
    var PA = P.multiply(this.data);
    for (var j = 0; j < n; ++j) {
      L[j*n+j] = 1;
      for (var i = 0; i < j+1; ++i) {
        var s = 0;
        for (var k = 0; k < i; ++k)
          s += U[k*n+j] * L[i*n+k]
        U[i*n+j] = PA[i*n+j] - s
      }
      for (var i = j; i < n; ++i) {
        var s = 0;
        for (var k = 0; k < i; ++k)
          s += U[k*n+j] * L[i*n+k]
        L[i*n+j] = (PA[i*n+j] - s) / U[j*n+j];
      }
    }
    return {
      L: new Array2d(L, n, n, true), 
      U: new Array2d(U, n, n, true), 
      P: P
    };
  }

  // A = LL^T
  cholInplace() {
    this.checkSquare();
    var A = this.data;
    var m = this.rows;
    var n = this.cols;
    var i, j, k;
    var s = 0.0;
    for (i = 0; i < n; ++i) {
      for (j = 0; j < (i + 1); ++j) {
        s = 0.0;
        for (k = 0; k < j; ++k)
          s += A[i * n + k] * A[j * n + k];
        if (i != j) A[j * n + i] = 0;
        if (i == j && A[i * n + i] - s < 0) throw "matrix not positive definite";
        A[i * n + j] = (i == j) ? Math.sqrt(A[i * n + i] - s) : ((A[i * n + j] - s) / A[j * n + j]);
      }
    }
    return Array2d(A, m, n, true);
  }

  // A = LL^T
  chol() {
    return this.copy().cholInplace();
  }

  // solve Lx = b using forward substitution, updates b
  forwardSubInplace(b) {
    var L = this.data;
    var m = this.rows, n = this.cols;
    for (var i = 0; i < n; ++i) {
      var s = 0.0
      for (var j = 0; j < i; ++j)
        s += L[i * n + j] * b.data[j];
      b.data[i] = (b.data[i] - s) / L[i * n + i];
    }
    return b;    
  }

  forwardSub(b) {
    return this.forwardSubInplace(b.copy());
  }

  // solve Ux = b using backward substitution, updates b
  backwardSubInplace(b, transpose) {
    var U = this.data;
    var m = this.rows, n = this.cols;
    for (var i = n - 1; i >= 0; --i) {
      var s = 0.0;
      for (var j = i + 1; j < n; ++j)
        s += (transpose ? U[j * n + i] : U[i * n + j]) * b.data[j];
      b.data[i] = (b.data[i] - s) / U[i * n + i];
    }
    return b;
  }

  backwardSub(b) {
    return this.backwardSubInplace(b.copy());
  }

  solve(b, method) {
    method = method || 'lu';
    if (method == 'lu') {
      let res = this.lu();
      let P = res.P, L = res.L, U = res.U;
      return U.backwardSubInplace(L.forwardSubInplace(P.multiply(b)));
    }
    if (method == 'chol') {
      var L = this.cholInplace();
      return L.backwardSubInplace(L.forwardSubInplace(b.copy()), true);
    }
    if (method == 'qr') {
      let res = this.qr();
      return res.R.backwardSubInplace(res.Q.transpose().multiply(b));
    }
  }

  inverse(b, method) {
    method = method || 'lu';
    if (method == 'lu') {
      var res = this.lu();
      var P = res.P, L = res.L, U = res.U;
      var inverse = Array2d.zeros(this.rows, this.cols);
      var eye = Array2d.eye(this.rows, this.cols);
      for (var j = 0; j < this.cols; ++j) {
        inverse.setCol(j, U.backwardSubInplace(L.forwardSubInplace(P.multiply(eye.col(j)))));
      }
      return inverse;
    }
    if (method == 'chol') {
      var L = this.chol();
      var inverse = Array2d.zeros(this.rows, this.cols);
      var eye = Array2d.eye(this.rows, this.cols);
      for (var j = 0; j < this.cols; ++j) {
        inverse.setCol(j, L.backwardSubInplace(L.forwardSubInplace(eye.col(j)), true));
      }
      return inverse;
    }
  }

  eig(options) {
    this.checkSquare();

    if (arguments.length < 1)
      options = {};

    var maxIter = options.maxIter || 100;
    var tolerance = options.tolerance || 1e-5;
    var verbose = options.verbose || false;

    var n = this.rows;
    var D = this.copy().data;
    var V = Array2d.eye(n, n);

    var iter, maxOffDiag, p, q;
    for (iter = 0; iter < maxIter; ++iter) {

      // find max off diagonal term at (p, q)
      maxOffDiag = 0;
      for (var i = 0; i < n - 1; ++i) {
        for (var j = i + 1; j < n; ++j) {
          if (Math.abs(D[i * n + j]) > maxOffDiag) {
            maxOffDiag = Math.abs(D[i * n + j]);
            p = i; q = j;
          }
        }
      }

      if (maxOffDiag < tolerance)
        break;

      // Rotates matrix D through theta in pq-plane to set D[p][q] = 0
      // Rotation stored in matrix V whose columns are eigenvectors of D
      // d = cot 2 * theta, t = tan theta, c = cos theta, s = sin theta
      var d = (D[p * n + p] - D[q * n + q]) / (2.0 * D[p * n + q]);
      var t = Math.sign(d) / (Math.abs(d) + Math.sqrt(d * d + 1));
      var c = 1.0 / Math.sqrt(t * t + 1);
      var s = t * c;
      D[p * n + p] += t * D[p * n + q];
      D[q * n + q] -= t * D[p * n + q];
      D[p * n + q] =      D[q * n + p] = 0.0;
      for (var k = 0; k < n; k++) {  // Transform D
        if (k != p && k != q) {
          var akp =  c * D[k * n + p] + s * D[k * n + q];
          var akq = -s * D[k * n + p] + c * D[k * n + q];
          D[k * n + p] = akp;
          D[p * n + k] = akp;
          D[k * n + q] = akq;
          D[q * n + k] = akq;
        }
      }
      for (var k = 0; k < n; k++) {  // Store V
        var rkp =  c * V.data[k * n + p] + s * V.data[k * n + q];
        var rkq = -s * V.data[k * n + p] + c * V.data[k * n + q];
        V.data[k * n + p] = rkp;
        V.data[k * n + q] = rkq;
      }
    }

    if (verbose && iter == maxIter) {
      console.log('Reached maxIter: ', maxOffDiag, ' > ', tolerance);
    }

    return { V: V, D: new Array2d(D, n, n, true) };

  }

  qr() {
    var makeHouseholder = function(a) {
      var v = a.scale(1 / (a.data[0] + Math.sign(a.data[0]) * a.norm()));
      v.data[0] = 1;
      var H = Array2d.eye(a.length);
      H.decrement(Array2d.outer(v, v).scale(2 / v.dot(v)));
      return H;
    };
    var A = this.copy();
    var m = A.rows;
    var n = A.cols;
    var Q = Array2d.eye(m);
    var upper = n - ((m == n) ? 1 : 0);
    for (var i = 0; i < upper; i++) {
      var a = A.getBlock(i, i, m - i, 1);
      var H = Array2d.eye(m);
      H.setBlock(i, i, makeHouseholder(a));
      Q = Q.multiply(H);
      A = H.multiply(A);
    }
    return {Q:Q, R:A};
  }

  toString(precision) {
    precision = precision || 4;
    let str = '';
    for (let i = 0; i < this.rows; i++) {
      str += (i == 0) ? '[[ ' : ' [ ';
      str += this.data[i * this.cols + 0].toPrecision(precision);
      for (let j = 1; j < this.cols; j++)
        str += ', ' + this.data[i * this.cols + j].toPrecision(precision);
      str += (i == this.rows - 1) ? ' ]]' : ' ],\n';
    }
    return str;    
  }

}

