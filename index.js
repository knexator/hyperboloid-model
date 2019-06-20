"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});

const mats = require('gl-matrix')

mats.mat3.getRow = function(out, mat, i) {
    return mats.vec3.set(out, mat[0+i], mat[3+i], mat[6+i]);
}

mats.mat3.setRow = function(out, vec, i) {
    out[0+i] = vec[0];
    out[3+i] = vec[1];
    out[6+i] = vec[2];
    return out;
}

exports.origin = mats.vec3.fromValues(0,0,1);
exports.metric = mats.mat3.create();
mats.mat3.fromScaling(exports.metric, mats.vec3.fromValues(1,1,-1));

exports.add_z = function(x, y, type=-1) {
    return mats.vec3.fromValues(x, y, Math.sqrt(Math.abs(type - x*x - y*y)));
}

exports.dot = function(pointA, pointB) {
    return pointA[0]*pointB[0] + pointA[1]*pointB[1] - pointA[2]*pointB[2];
    /*var p = mats.vec3.create();
    mats.vec3.transformMat3(p, pointA, exports.metric);
    var q = mats.vec3.create();
    mats.vec3.transformMat3(q, pointB, exports.metric);
    return mats.vec3.dot(p,q);*/
}

exports.norm = function(point) {
    return Math.sqrt(Math.abs(exports.dot(point, point)));
}

exports.normalize = function(out, point) {
    var norm = exports.norm(point);
    if (norm == 0) {
        if (point[2] != 0) {
            out[0] = point[0]/point[2];
            out[1] = point[1]/point[2];
            out[2] = 1; //point[2]/point[2];
            return out;
        } else {
            out[0] = 0;
            out[1] = 0;
            out[2] = 0;
            return out;
        }
    } else {
        out[0] = point[0]/norm;
        out[1] = point[1]/norm;
        out[2] = point[2]/norm;
        return out;
    }
}

exports.type = function(a) {
    // -1 is a point, 0 is an ideal point, 1 is a line
    return Math.sign(exports.dot(a, a));
}

exports.is_point = function(point) {
    return exports.dot(point, point) < 0;
}

exports.is_line = function(point) {
    return exports.dot(point, point) > 0;
}

exports.is_ideal = function(point) {
    return exports.dot(point, point) == 0;
}

exports.line_joining_points = function(out, pointA, pointB) {
    mats.vec3.cross(out, pointA, pointB);
    out[2] = -out[2];
    return exports.normalize(out, out);
}

exports.ideals_at_infinity = function(outA, outB, line) {
  /*# the target lines must obey
  # -   w*w +   x*x +   y*y = 0, (belong to the cone) and
  # - v.w*w + v.x*x + v.y*y = 0  (belong to the dual of v)
  
  # taking w = 1, we get x*x + y*y = 1, and
  # - v.w + v.x*x + v.y*y = 0 ->
  # x = (v.w - v.y*y)/v.x (if v.x is not 0; if it is, we do it with v.y)
  # so then
  # ((v.w - v.y*y)/v.x)**2 + y**2 = 1  gives us y,
  # v.w**2 - 2*v.w*v.y*y + (v.x*v.x+v.y*v.y)*y*y = v.x*v.x*/
    if (exports.type(line) == -1) {
        throw new Error("line is a point, in ideals_at_infinity");
    }
    
    let line_n = mats.vec3.fromValues(line[0], line[1], -line[2]);
    exports.normalize(line_n, line_n);
    if (line_n[0] != 0) {
        let a = line_n[0]*line_n[0]+line_n[1]*line_n[1];
        let b = 2*line_n[1]*line_n[2];
        let c = line_n[2]*line_n[2] - line_n[0]*line_n[0];
        let disc = b*b - 4*a*c;
        if (disc < 0) {
            throw new Error("negative discriminant at ideals_at_infinity");
        } else {
            let order = Math.sign(line_n[0]);
            let y1 = (-b + order*Math.sqrt(disc)) / (2*a);
            let y2 = (-b - order*Math.sqrt(disc)) / (2*a);
            let x1 = (-y1 * line_n[1] - line_n[2]) / line_n[0];
            let x2 = (-y2 * line_n[1] - line_n[2]) / line_n[0];            
            return [mats.vec3.set(outA, x1, y1, 1), mats.vec3.set(outB, x2, y2, 1)];
        }
    } else if (line_n[1] != 0) {
        let a = line_n[0]*line_n[0]+line_n[1]*line_n[1];
        let b = 2*line_n[0]*line_n[2];
        let c = line_n[2]*line_n[2] - line_n[1]*line_n[1];
        let disc = b*b - 4*a*c;
        if (disc < 0) {
            throw new Error("negative discriminant at ideals_at_infinity");
        } else {
            let order = Math.sign(line_n[1]);
            let x1 = (-b + order*Math.sqrt(disc)) / (2*a);
            let x2 = (-b - order*Math.sqrt(disc)) / (2*a);
            let y1 = (-x1 * line_n[0] - line_n[2]) / line_n[1];
            let y2 = (-x2 * line_n[0] - line_n[2]) / line_n[1];            
            return [mats.vec3.set(outA, x1, y1, 1), mats.vec3.set(outB, x2, y2, 1)];
        }
    } else {
        throw new Error("Unexpected value, in ideals_at_infinity");
    }
}

exports.dist = function(a, b) {
    var dot = exports.dot(a, b);
    var a_type = exports.type(a);
    var b_type = exports.type(b);
    
    if (a_type == -1 && b_type == -1) { // both points
        // Distance between two points
        return Math.acosh(-dot);
    }
    if (a_type == -1 && b_type == 1) { // point and line
        return Math.asinh(dot);
    }
    if (a_type == 1 && b_type == 1) { // both lines
        if (Math.abs(dot) > 1) {
            // ultraparallel, return distance between them
            return Math.acosh(Math.abs(dot));
        } else if (Math.abs(dot) == 1) {
            // parallel, they meet at infinity
            return 0;            
        } else {
            // both lines intersect, return angle between them
            return Math.acos(dot);
        }
    }
    // TODO: handle cases for ideal points
}

exports.hlerp = function(out, pointA, pointB, t) {
    var dist = exports.dist(pointA, pointB);
    var a = Math.sinh((1-t) * dist);
    var b = Math.sinh(t * dist);
    out[0] = a * pointA[0] + b * pointB[0];
    out[1] = a * pointA[1] + b * pointB[1];
    out[2] = a * pointA[2] + b * pointB[2];
    exports.normalize(out, out);
    return out;
}

exports.rot = function(out, ang) {
    let c = Math.cos(ang);
    let s = Math.sin(ang);
    return mats.mat3.set(out, c, s, 0, -s, c, 0, 0, 0, 1);
}

exports.xMove = function(out, dist) {
    let c = Math.cosh(dist);
    let s = Math.sinh(dist);
    return mats.mat3.set(out, c, 0, s, 0, 1, 0, s, 0, c);
}

exports.translation_along_line = function(out, line, distance) {
    if (distance == 0) {
        return mats.mat3.identity(out);
    }
    let lineA = mats.vec3.create();
    exports.normalize(lineA, line);
    let lineB = mats.vec3.create();
    let lineC = mats.vec3.create();
    exports.ideals_at_infinity(lineB, lineC, line);
    
    let c = Math.cosh(distance);
    let s = Math.sinh(distance);
    let M = mats.mat3.fromValues(1, 0, 0,  0, c-s, 0,  0, 0, c+s);
    let S = mats.mat3.fromValues(...lineA, ...lineB, ...lineC);
    let R = mats.mat3.create();
    mats.mat3.invert(R, S);
    mats.mat3.mul(M, M, R);
    return mats.mat3.mul(out, S, M);
}

exports.translation_between_points = function(out, pointA, pointB) {
    if (mats.vec3.equals(pointA, pointB)) {
        return mats.mat3.identity(out);
    }
    let line = mats.vec3.create();
    exports.line_joining_points(line, pointA, pointB);
    let distance = exports.dist(pointA, pointB);
    return exports.translation_along_line(out, line, distance);
}

exports.gram_schmidt = function(out, mat) {
    mats.mat3.copy(out, mat);
    for (let i=0; i<3; i++) {
        let rowI = mats.mat3.getRow(mats.vec3.create(), out, i);
        let rowNorm = exports.norm(rowI);
        if (rowNorm == 0) {
            throw new Error("velocity row");
        }
        //console.log(mats.vec3.scale(rowI, rowI, 1.0/rowNorm));
        mats.vec3.scale(rowI, rowI, 1.0/rowNorm);
        mats.mat3.setRow(out, rowI, i);
        for (let j=i+1; j<3; j++) {
            let rowJ = mats.mat3.getRow(mats.vec3.create(), out, j);
            let component = exports.dot(rowI, rowJ);
            mats.vec3.scaleAndAdd(rowJ, rowJ, rowI, -component);
            mats.mat3.setRow(out, rowJ, j);
        }
    }
    return out;
}