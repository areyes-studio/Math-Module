class Vector3 {
	/**
	 * Creates an instance of Vector3.
	 * @param {number} x
	 * @param {number} y
	 * @param {number} z
	 * @memberof Vector3
	 */
	constructor(x, y, z) {
		Object.defineProperty(this, "x", { get: () => x });
		Object.defineProperty(this, "y", { get: () => y });
		Object.defineProperty(this, "z", { get: () => z });
	}

	get x() { return 0; }
	get y() { return 0; }
	get z() { return 0; }

	static get zero() { return new Vector3(0, 0, 0); }
	static get right() { return new Vector3(1, 0, 0); }
	static get up() { return new Vector3(0, 1, 0); }
	static get forward() { return new Vector3(0, 0, 1); }

	static add(/** @type {Vector3} */ vec1, /** @type {Vector3} */ vec2) {
		return new Vector3(
			vec1.x + vec2.x,
			vec1.y + vec2.y,
			vec1.z + vec2.z
		);
	};

	static sub(/** @type {Vector3} */ vec1, /** @type {Vector3} */ vec2) {
		return new Vector3(
			vec1.x - vec2.x,
			vec1.y - vec2.y,
			vec1.z - vec2.z
		);
	};

	static mul(/** @type {Vector3} */ vec1, /** @type {number} */ scalar) {
		return new Vector3(
			vec1.x * scalar,
			vec1.y * scalar,
			vec1.z * scalar
		);
	};

	static div(/** @type {Vector3} */ vec1, /** @type {number} */ scalar) {
		return new Vector3(
			vec1.x / scalar,
			vec1.y / scalar,
			vec1.z / scalar
		);
	};

	/**
	 * Returns dot product of two vectors
	 *
	 * @static
	 * @param {Vector3} vec1
	 * @param {Vector3} vec2
	 * @returns {number}
	 * @memberof Vector3
	 */
	static dot(vec1, vec2) {
		return vec1.x * vec2.x + vec1.y * vec2.y + vec1.z * vec2.z;
	};

	/**
	 * Returns cross product of two vectors
	 *
	 * @static
	 * @param {Vector3} vec1
	 * @param {Vector3} vec2
	 * @returns
	 * @memberof Vector3
	 */
	static cross(vec1, vec2) {
		return new Vector3(
			vec1.y * vec2.z - vec1.z * vec2.y,
			vec1.z * vec2.x - vec1.x * vec2.z,
			vec1.x * vec2.y - vec1.y * vec2.x
		);
	}

	/**
	 * Returns an angle between two vectors in radians
	 *
	 * @static
	 * @param {Vector3} vec1
	 * @param {Vector3} vec2
	 * @returns
	 * @memberof Vector3
	 */
	static angleBetween(vec1, vec2) {
		return Math.acos(Vector3.dot(vec1.normalize(), vec2.normalize()) / (vec1.normalize().magnitude * vec2.normalize().magnitude));
	}

	/**
	 * Returns signed angle between two vectors on specified axis
	 *
	 * @static
	 * @param {Vector3} vec1
	 * @param {Vector3} vec2
	 * @param {Vector3} axis
	 * @memberof Vector3
	 */
	static angleBetweenSigned(vec1, vec2, axis) {
		let unsignedAngle = Vector3.angleBetween(vec1, vec2);

		let cross_x = vec1.y * vec2.z - vec1.z * vec2.y;
		let cross_y = vec1.z * vec2.x - vec1.x * vec2.z;
		let cross_z = vec1.x * vec2.y - vec1.y * vec2.x;
		let sign = Math.sign(axis.x * cross_x + axis.y * cross_y + axis.z * cross_z);
		return unsignedAngle * sign;
	}

	/**
	 * Returns the number which is the distance between two vectors
	 * @param {Vector3} vector1
	 * @param {Vector3} vector2 
	 */
	static distance(vector1, vector2) {
		let subVec = new Vector3(
			vector1.x - vector2.x,
			vector1.y - vector2.y,
			vector1.z - vector2.z
		);

		return Math.sqrt(subVec.x ** 2 + subVec.y ** 2 + subVec.z ** 2);
	}

	/**
	 * Rotates vector by quaternion
	 *
	 * @static
	 * @param {*} v
	 * @param {*} q
	 * @returns
	 * @memberof Vector3
	 */
	static rotate(/** @type {Vector3} */ v, /** @type {{ w: number, x: number, y: number, z: number }} */ q) {
		let u = new Vector3(q.x, q.y, q.z);
		let s = q.w;

		return u.mul(2.0 * Vector3.dot(u, v))
			.add(v.mul(s * s - Vector3.dot(u, u)))
			.add(Vector3.cross(u, v).mul(2.0 * s));
	}

	/**
	 * Projects a vector onto a plane defined by normal
	 *
	 * @static
	 * @param {Vector3} vector
	 * @param {Vector3} planeNormal
	 * @memberof Vector3
	 */
	static projectOnPlane(vector, planeNormal) {
		let sqrMag = Vector3.dot(planeNormal, planeNormal);
		if (sqrMag < Number.EPSILON)
			return vector;
		else {
			let dot = Vector3.dot(vector, planeNormal);
			return new Vector3(
				vector.x - planeNormal.x * dot / sqrMag,
				vector.y - planeNormal.y * dot / sqrMag,
				vector.z - planeNormal.z * dot / sqrMag
			);
		}
	}

	/**
	 * Linearly interpolates between two vectors
	 *
	 * @param {Vector3} vec1
	 * @param {Vector3} vec2
	 * @param {number} t Interpolation factor - should be in range [0; 1]
	 * @return {Vector3} 
	 */
	static lerp(vec1, vec2, t) {
		t = AMath$1.clamp(t, 0, 1);
		return new Vector3(
			vec1.x + (vec2.x - vec1.x) * t,
			vec1.y + (vec2.y - vec1.y) * t,
			vec1.z + (vec2.z - vec1.z) * t
		);
	}

	add(/** @type {Vector3} */ vec) {
		return Vector3.add(this, vec);
	};

	sub(/** @type {Vector3} */ vec) {
		return Vector3.sub(this, vec);
	};

	/**
	 *
	 *
	 * @param {number} scalar
	 * @memberof Vector3
	 */
	mul(scalar) {
		return Vector3.mul(this, scalar);
	}

	/**
	 *
	 *
	 * @param {number} scalar
	 * @memberof Vector3
	 */
	div(scalar) {
		return Vector3.div(this, scalar);
	}

	/**
	 * Returns dot product of two vectors
	 *
	 * @static
	 * @param {Vector3} vec
	 * @returns {number}
	 * @memberof Vector3
	 */
	dot(vec) {
		return Vector3.dot(this, vec)
	};

	/**
	 * Returns cross product of two vectors
	 *
	 * @static
	 * @param {Vector3} vec
	 * @returns
	 * @memberof Vector3
	 */
	cross(vec) {
		return Vector3.cross(this, vec);
	}

	/**
	 * Returns an angle between two vectors in radians
	 *
	 * @static
	 * @param {Vector3} vec
	 * @returns
	 * @memberof Vector3
	 */
	angleBetween(vec) {
		return Vector3.angleBetween(this, vec);
	}

	/**
	 * Returns signed angle between two vectors on specified axis
	 *
	 * @static
	 * @param {Vector3} vector
	 * @param {Vector3} axis
	 * @memberof Vector3
	 */
	angleBetweenSigned(vector, axis) {
		return Vector3.angleBetweenSigned(this, vector, axis);
	}

	/**
	 * Returns the number which is the distance between two vectors
	 * @param {Vector3} vector 
	 */
	distance(vector) {
		return Vector3.distance(this, vector);
	}

	rotate(/** @type {{ w: number, x: number, y: number, z: number }} */ q) {
		return Vector3.rotate(this, q);
	}

	/**
	 * Returns the vector, the madnitude of which is limited by specified one
	 *
	 * @param {number} maxLength
	 * @returns
	 */
	clampMagnitude(maxLength) {
		let sqrmag = this.magnitudeSquared;
		if (sqrmag > maxLength * maxLength) {
			let mag = Math.sqrt(sqrmag);
			let normalized_x = this.x / mag;
			let normalized_y = this.y / mag;
			let normalized_z = this.z / mag;
			
			return new Vector3(
				normalized_x * maxLength,
				normalized_y * maxLength,
				normalized_z * maxLength
			);
		}
		return this;
	}

	/**
	 * Projects a vector onto a plane defined by normal
	 *
	 * @param {Vector3} planeNormal
	 * @memberof Vector3
	 */
	projectOnPlane(planeNormal) {
		return Vector3.projectOnPlane(this, planeNormal);
	}

	/**
	 * Linearly interpolates between two vectors
	 *
	 * @param {Vector3} vector
	 * @param {number} t Interpolation factor - should be in range [0; 1]
	 * @return {Vector3} 
	 */
	lerp(vector, t) {
		return Vector3.lerp(this, vector, t)
	}

	neg() {
		return this.mul(-1);
	}

	get magnitude() {
		return Math.sqrt(this.x * this.x + this.y * this.y + this.z * this.z);
	}

	get magnitudeSquared() {
		return this.x * this.x + this.y * this.y + this.z * this.z;
	}

	normalize() {
		return this.div(this.magnitude);
	}

	/**
	 * @static
	 * @param {PointSignal} point
	 * @returns
	 * @memberof Vector3
	 */
	static fromPointSignal(point) {
		return new Vector3(
			point.x.pinLastValue(),
			point.y.pinLastValue(),
			point.z.pinLastValue()
		);
	}

	/**
	 *
	 *
	 * @static
	 * @param {number[]} array
	 * @returns
	 * @memberof Vector3
	 */
	static fromArray(array) {
		return new Vector3(
			array[0],
			array[1],
			array[2]
		)
	}

	toPointSignal() {
		const Reactive = require("Reactive");
		return Reactive.point(this.x, this.y, this.z);
	}

	toVectorSignal() {
		const Reactive = require("Reactive");
		return Reactive.vector(this.x, this.y, this.z);
	}

	toString(digits = 2) {
		return "(" + this.x.toFixed(digits) + "; " + this.y.toFixed(digits) + "; " + this.z.toFixed(digits) + ")";
	}

	/**
	 * @return {[number, number, number]} 
	 * @memberof Vector3
	 */
	toArray() {
		return [this.x, this.y, this.z];
	}
}

class Quaternion {
    /**
     * Creates an instance of Vector3.
     * @param {number} w
     * @param {number} x
     * @param {number} y
     * @param {number} z
     * @memberof Quaternion
     */
    constructor(w, x, y, z) {
        if (Array.isArray(w)) {
            [w, x, y, z] = w;
        } else if (Array.isArray(x)) {
            [x, y, z] = x;
        } else if (Array.isArray(y)) {
            [y, z] = y;
        }

        Object.defineProperty(this, "w", { get: () => w ? w : 1});
        Object.defineProperty(this, "x", { get: () => x ? x : 0});
        Object.defineProperty(this, "y", { get: () => y ? y : 0});
        Object.defineProperty(this, "z", { get: () => z ? z : 0});
    }

    get w() { return 1; }
    get x() { return 0; }
    get y() { return 0; }
    get z() { return 0; }

    static get zero() { return new Quaternion(0, 0, 0, 0); }
    static get one() { return new Quaternion(1, 0, 0, 0); }
    static get i() { return new Quaternion(0, 1, 0, 0); }
    static get j() { return new Quaternion(0, 0, 1, 0); }
    static get k() { return new Quaternion(0, 0, 0, 1); }
    static get epsilon() { return 1e-16 }

    /**
     * @static
     * @param {number[]} axis
     * @param {number} angle
     * @return {Quaternion} 
     * @memberof Quaternion
     */
    static fromAxisAngle(axis, angle) {

        // Q = [cos(angle / 2), v * sin(angle / 2)]

        var halfAngle = angle * 0.5;

        var a = axis[0];
        var b = axis[1];
        var c = axis[2];

        var sin_2 = Math.sin(halfAngle);
        var cos_2 = Math.cos(halfAngle);

        var sin_norm = sin_2 / Math.sqrt(a * a + b * b + c * c);

        return new Quaternion(cos_2, a * sin_norm, b * sin_norm, c * sin_norm);
    }

    /**
     *
     * @static
     * @param {Array<number>} u
     * @param {Array<number>} v
     * @return {*} 
     * @memberof Quaternion
     */
    static fromBetweenVectors(u, v) {
        var ux = u[0];
        var uy = u[1];
        var uz = u[2];

        var vx = v[0];
        var vy = v[1];
        var vz = v[2];

        var dot = ux * vx + uy * vy + uz * vz;

        // Parallel check (TODO must be normalized)
        if (dot >= 1 - Quaternion.epsilon) ;

        // Close to PI @TODO
        //if (1 + dot <= Quaternion['EPSILON']) {
        // return Quaternion.fromAxisAngle(Math.abs(ux) > Math.abs(uz) ? [-uy,  ux, 0] : [0, -uz,  uy], 0) -> OR
        // return Quaternion.fromAxisAngle(Math.abs(ux) > Math.abs(uz) ? [ uy, -ux, 0] : [0,  uz, -uy], 0)
        //}

        var wx = uy * vz - uz * vy;
        var wy = uz * vx - ux * vz;
        var wz = ux * vy - uy * vx;

        return new Quaternion(
            dot + Math.sqrt(dot * dot + wx * wx + wy * wy + wz * wz),
            wx,
            wy,
            wz
        ).normalize();
    }

    /**
     * @static
     * @return {*} 
     * @memberof Quaternion
     */
    static random() {

        var u1 = Math.random();
        var u2 = Math.random();
        var u3 = Math.random();

        var s = Math.sqrt(1 - u1);
        var t = Math.sqrt(u1);

        return new Quaternion(
            t * Math.cos(2 * Math.PI * u3),
            s * Math.sin(2 * Math.PI * u2),
            s * Math.cos(2 * Math.PI * u2),
            t * Math.sin(2 * Math.PI * u3)
        );
    }

    /**
     * @static
     * @param {number} phi
     * @param {number} theta
     * @param {number} psi
     * @param {string=} order
     * @memberof Quaternion
     */
    static fromEuler(phi, theta, psi, order) {
        var _x = phi * 0.5;
        var _y = theta * 0.5;
        var _z = psi * 0.5;

        var cX = Math.cos(_x);
        var cY = Math.cos(_y);
        var cZ = Math.cos(_z);

        var sX = Math.sin(_x);
        var sY = Math.sin(_y);
        var sZ = Math.sin(_z);

        if (order === undefined || order === 'ZXY') {
            // axisAngle([0, 0, 1], x) * axisAngle([1, 0, 0], y) * axisAngle([0, 1, 0], z)
            return new Quaternion(
                cX * cY * cZ - sX * sY * sZ,
                cX * cZ * sY - cY * sX * sZ,
                cX * cY * sZ + cZ * sX * sY,
                cY * cZ * sX + cX * sY * sZ);
        }

        if (order === 'XYZ' || order === 'RPY') {
            // axisAngle([1, 0, 0], x) * axisAngle([0, 1, 0], y) * axisAngle([0, 0, 1], z)
            return new Quaternion(
                cX * cY * cZ - sX * sY * sZ,
                cY * cZ * sX + cX * sY * sZ,
                cX * cZ * sY - cY * sX * sZ,
                cX * cY * sZ + cZ * sX * sY);
        }

        if (order === 'YXZ') {
            // axisAngle([0, 1, 0], x) * axisAngle([1, 0, 0], y) * axisAngle([0, 0, 1], z)
            return new Quaternion(
                cX * cY * cZ + sX * sY * sZ,
                cX * cZ * sY + cY * sX * sZ,
                cY * cZ * sX - cX * sY * sZ,
                cX * cY * sZ - cZ * sX * sY);
        }

        if (order === 'ZYX' || order === 'YPR') {
            // axisAngle([0, 0, 1], x) * axisAngle([0, 1, 0], y) * axisAngle([1, 0, 0], z)
            return new Quaternion(
                cX * cY * cZ + sX * sY * sZ,
                cX * cY * sZ - cZ * sX * sY,
                cX * cZ * sY + cY * sX * sZ,
                cY * cZ * sX - cX * sY * sZ);
        }

        if (order === 'YZX') {
            // axisAngle([0, 1, 0], x) * axisAngle([0, 0, 1], y) * axisAngle([1, 0, 0], z)
            return new Quaternion(
                cX * cY * cZ - sX * sY * sZ,
                cX * cY * sZ + cZ * sX * sY,
                cY * cZ * sX + cX * sY * sZ,
                cX * cZ * sY - cY * sX * sZ);
        }

        if (order === 'XZY') {
            // axisAngle([1, 0, 0], x) * axisAngle([0, 0, 1], z) * axisAngle([0, 1, 0], y)
            return new Quaternion(
                cX * cY * cZ + sX * sY * sZ,
                cY * cZ * sX - cX * sY * sZ,
                cX * cY * sZ - cZ * sX * sY,
                cX * cZ * sY + cY * sX * sZ);
        }
        return null;
    }

    /**
     * @static
     * @param {Array<*>} matrix
     * @memberof Quaternion
     */
    static fromMatrix(matrix) {
        if (matrix.length === 9) {

            var m00 = matrix[0];
            var m01 = matrix[1];
            var m02 = matrix[2];

            var m10 = matrix[3];
            var m11 = matrix[4];
            var m12 = matrix[5];

            var m20 = matrix[6];
            var m21 = matrix[7];
            var m22 = matrix[8];

        } else {
            var m00 = matrix[0][0];
            var m01 = matrix[0][1];
            var m02 = matrix[0][2];

            var m10 = matrix[1][0];
            var m11 = matrix[1][1];
            var m12 = matrix[1][2];

            var m20 = matrix[2][0];
            var m21 = matrix[2][1];
            var m22 = matrix[2][2];
        }

        var tr = m00 + m11 + m22;

        if (tr > 0) {
            var S = Math.sqrt(tr + 1.0) * 2; // S=4*qw

            return new Quaternion(
                0.25 * S,
                (m21 - m12) / S,
                (m02 - m20) / S,
                (m10 - m01) / S);
        } else if ((m00 > m11) && (m00 > m22)) {
            var S = Math.sqrt(1.0 + m00 - m11 - m22) * 2; // S=4*qx

            return new Quaternion(
                (m21 - m12) / S,
                0.25 * S,
                (m01 + m10) / S,
                (m02 + m20) / S);
        } else if (m11 > m22) {
            var S = Math.sqrt(1.0 + m11 - m00 - m22) * 2; // S=4*qy

            return new Quaternion(
                (m02 - m20) / S,
                (m01 + m10) / S,
                0.25 * S,
                (m12 + m21) / S);
        } else {
            var S = Math.sqrt(1.0 + m22 - m00 - m11) * 2; // S=4*qz

            return new Quaternion(
                (m10 - m01) / S,
                (m02 + m20) / S,
                (m12 + m21) / S,
                0.25 * S);
        }
    }

    /**
     * Adds two quaternions Q1 and Q2
     *
     * @static
     * @param {Quaternion} quat1
     * @param {Quaternion} quat2
     * @return {*} 
     * @memberof Quaternion
     */
    static add(quat1, quat2) {
        return new Quaternion(
            quat1.w + quat2.w,
            quat1.x + quat2.x,
            quat1.y + quat2.y,
            quat1.z + quat2.z
        );
    };

    /**
     * Subtracts two quaternions Q1 and Q2
     *
     * @static
     * @param {Quaternion} quat1
     * @param {Quaternion} quat2
     * @return {*} 
     * @memberof Quaternion
     */
    static sub(/** @type {Quaternion} */ quat1, /** @type {Quaternion} */ quat2) {
        return new Quaternion(
            quat1.w - quat2.w,
            quat1.x - quat2.x,
            quat1.y - quat2.y,
            quat1.z - quat2.z
        );
    };

    /**
     * Calculates the additive inverse, or simply it negates the quaternion
     *
     * @static
     * @param {Quaternion} quat
     * @return {*} 
     * @memberof Quaternion
     */
    static neg(quat) {
        return new Quaternion(
            -quat.w,
            -quat.x,
            -quat.y,
            -quat.z
        );
    };

    /**
     * Calculates the length/modulus/magnitude or the norm of a quaternion
     *
     * @static
     * @param {Quaternion} quat
     * @return {*} 
     * @memberof Quaternion
     */
    static norm(quat) {
        return Math.sqrt(quat.w * quat.w + quat.x * quat.x + quat.y * quat.y + quat.z * quat.z)
    };

    /**
  * Calculates the squared length/modulus/magnitude or the norm of a quaternion
     * @static
     * @param {Quaternion} quat
     * @return {*} 
     * @memberof Quaternion
     */
    static normSq(quat) {
        return quat.w * quat.w + quat.x * quat.x + quat.y * quat.y + quat.z * quat.z
    };

    /**
     * Normalizes the quaternion to have |Q| = 1 as long as the norm is not zero
     * Alternative names are the signum, unit or versor
     * @static
     * @param {Quaternion} quat
     * @return {*} 
     * @memberof Quaternion
     */
    static normalize(quat) {
        var norm = Math.sqrt(quat.w * quat.w + quat.x * quat.x + quat.y * quat.y + quat.z * quat.z);

        if (norm < Quaternion.epsilon) {
            return Quaternion.zero;
        }

        norm = 1 / norm;

        return new Quaternion(quat.w * norm, quat.x * norm, quat.y * norm, quat.z * norm);
    };

    /**
     * Calculates the Hamilton product of two quaternions
     * Leaving out the imaginary part results in just scaling the quat
     * @static
     * @param {*} quat1
     * @param {*} quat2
     * @return {*} 
     * @memberof Quaternion
     */
    static mul(/** @type {Quaternion} */ quat1, /** @type {Quaternion} */ quat2) {
        var w1 = quat1.w;
        var x1 = quat1.x;
        var y1 = quat1.y;
        var z1 = quat1.z;

        var w2 = quat2.w;
        var x2 = quat2.x;
        var y2 = quat2.y;
        var z2 = quat2.z;

        return new Quaternion(
            w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2,
            w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2,
            w1 * y2 + y1 * w2 + z1 * x2 - x1 * z2,
            w1 * z2 + z1 * w2 + x1 * y2 - y1 * x2
        );
    };

    /**
     * Scales a quaternion by a scalar, faster than using multiplication
     *
     * @static
     * @param {Quaternion} quat
     * @param {number} s
     * @return {*} 
     * @memberof Quaternion
     */
    static scale(quat, s) {
        return new Quaternion(
            quat.w * s,
            quat.x * s,
            quat.y * s,
            quat.z * s
        );
    };

    /**
     * Calculates the dot product of two quaternions
     *
     * @static
     * @param {Quaternion} quat1
     * @param {Quaternion} quat2
     * @return {*} 
     * @memberof Quaternion
     */
    static dot(quat1, quat2) {
        return quat1.w * quat2.w + quat1.x * quat2.x + quat1.y * quat2.y + quat1.z * quat2.z
    };

    /**
     * Calculates the inverse of a quat for non-normalized quats such that
     * Q^-1 * Q = 1 and Q * Q^-1 = 1
     *
     * @static
     * @param {Quaternion} quat
     * @return {*} 
     * @memberof Quaternion
     */
    static inverse(quat) {
        var w = quat.w;
        var x = quat.x;
        var y = quat.y;
        var z = quat.z;

        var normSq = w * w + x * x + y * y + z * z;

        if (normSq === 0) {
            return Quaternion.zero; // TODO: Is the result zero or one when the norm is zero?
        }

        normSq = 1 / normSq;

        return new Quaternion(w * normSq, -x * normSq, -y * normSq, -z * normSq);
    }

    /**
     * Multiplies a quaternion with the inverse of a second quaternion
     * @static
     * @param {Quaternion} quat1
     * @param {Quaternion} quat2
     * @return {*} 
     * @memberof Quaternion
     */
    static div(quat1, quat2) {

        var w1 = quat1.w;
        var x1 = quat1.x;
        var y1 = quat1.y;
        var z1 = quat1.z;

        var w2 = quat2.w;
        var x2 = quat2.x;
        var y2 = quat2.y;
        var z2 = quat2.z;

        var normSq = w2 * w2 + x2 * x2 + y2 * y2 + z2 * z2;

        if (normSq === 0) {
            return Quaternion.zero; // TODO: Is the result zero or one when the norm is zero?
        }

        normSq = 1 / normSq;
        return new Quaternion(
            (w1 * w2 + x1 * x2 + y1 * y2 + z1 * z2) * normSq,
            (x1 * w2 - w1 * x2 - y1 * z2 + z1 * y2) * normSq,
            (y1 * w2 - w1 * y2 - z1 * x2 + x1 * z2) * normSq,
            (z1 * w2 - w1 * z2 - x1 * y2 + y1 * x2) * normSq
        );
    };

    /**
     * Calculates the conjugate of a quaternion
     * @static
     * @param {Quaternion} quat
     * @return {*} 
     * @memberof Quaternion
     */
    static conjugate(quat) {
        return new Quaternion(
            quat.w,
            -quat.x,
            -quat.y,
            -quat.z
        );
    };

    /**
     *
     * @static
     * @param {Quaternion} quat
     * @return {*} 
     * @memberof Quaternion
     */
    static exp(quat) {
        var w = quat.w;
        var x = quat.x;
        var y = quat.y;
        var z = quat.z;

        var vNorm = Math.sqrt(x * x + y * y + z * z);
        var wExp = Math.exp(w);
        var scale = wExp / vNorm * Math.sin(vNorm);

        if (vNorm === 0) {
            //return new Quaternion(wExp * Math.cos(vNorm), 0, 0, 0);
            return new Quaternion(wExp, 0, 0, 0);
        }

        return new Quaternion(
            wExp * Math.cos(vNorm),
            x * scale,
            y * scale,
            z * scale);
    };

    /**
     * Calculates the natural logarithm of the quaternion
     * @static
     * @param {Quaternion} quat
     * @return {*} 
     * @memberof Quaternion
     */
    static log(quat) {
        var w = quat.w;
        var x = quat.x;
        var y = quat.y;
        var z = quat.z;

        if (y === 0 && z === 0) {
            return new Quaternion(
                logHypot(w, x),
                Math.atan2(x, w), 0, 0);
        }

        var qNorm2 = x * x + y * y + z * z + w * w;
        var vNorm = Math.sqrt(x * x + y * y + z * z);

        var scale = Math.atan2(vNorm, w) / vNorm;

        return new Quaternion(
            Math.log(qNorm2) * 0.5,
            x * scale,
            y * scale,
            z * scale);
    }

    /**
     * Calculates the power of a quaternion raised to a real number or another quaternion
     * @static
     * @param {Quaternion} quat1
     * @param {Quaternion} quat2
     * @return {*} 
     * @memberof Quaternion
     */
    static pow(quat1, quat2) {

        var w1 = quat1.w;
        var x1 = quat1.x;
        var y1 = quat1.y;
        var z1 = quat1.z;

        var w2 = quat2.w;
        var x2 = quat2.x;
        var y2 = quat2.y;
        var z2 = quat2.z;

        if (y2 === 0 && z2 === 0) {

            if (w2 === 1 && x2 === 0) {
                return quat1;
            }

            if (w2 === 0 && x2 === 0) {
                return Quaternion.one;
            }

            // Check if we can operate in C
            // Borrowed from complex.js
            if (y1 === 0 && z1 === 0) {

                var a = w1;
                var b = x1;

                if (a === 0 && b === 0) {
                    return Quaternion.zero;
                }

                var arg = Math.atan2(b, a);
                var loh = logHypot(a, b);

                if (x2 === 0) {

                    if (b === 0 && a >= 0) {

                        return new Quaternion(Math.pow(a, w2), 0, 0, 0);

                    } else if (a === 0) {

                        switch (w2 % 4) {
                            case 0:
                                return new Quaternion(Math.pow(b, w2), 0, 0, 0);
                            case 1:
                                return new Quaternion(0, Math.pow(b, w2), 0, 0);
                            case 2:
                                return new Quaternion(-Math.pow(b, w2), 0, 0, 0);
                            case 3:
                                return new Quaternion(0, -Math.pow(b, w2), 0, 0);
                        }
                    }
                }

                a = Math.exp(w2 * loh - x2 * arg);
                b = x2 * loh + w2 * arg;
                return new Quaternion(
                    a * Math.cos(b),
                    a * Math.sin(b), 0, 0);
            }
        }

        // Normal quaternion behavior
        // q^p = e^ln(q^p) = e^(ln(q)*p)
        return quat1.log().mul(quat2).exp();
    };

    /**
     * Checks if two quats are the same
     * @static
     * @param {Quaternion} quat1
     * @param {Quaternion} quat2
     * @return {*} 
     * @memberof Quaternion
     */
    static equals(quat1, quat2) {

        let eps = Quaternion.epsilon;

        return Math.abs(quat2.w - quat1.w) < eps
            && Math.abs(quat2.x - quat1.x) < eps
            && Math.abs(quat2.y - quat1.y) < eps
            && Math.abs(quat2.z - quat1.z) < eps;
    };

    /**
     * Checks if all parts of a quaternion are finite
     * @static
     * @param {Quaternion} quat
     * @return {*} 
     * @memberof Quaternion
     */
    static isFinite(quat) {

        return isFinite(quat.w) && isFinite(quat.x) && isFinite(quat.y) && isFinite(quat.z);
    };

    /**
     * Checks if any of the parts of the quaternion is not a number
     * @static
     * @param {Quaternion} quat
     * @return {*} 
     * @memberof Quaternion
     */
    static isNaN(quat) {
        return isNaN(quat.w) || isNaN(quat.x) || isNaN(quat.y) || isNaN(quat.z);
    };

    /**
     * Gets the Quaternion as a well formatted string
     * @static
     * @param {Quaternion} quat
     * @return {*} 
     * @memberof Quaternion
     */
    static toString(quat) {
        var w = quat.w;
        var x = quat.x;
        var y = quat.y;
        var z = quat.z;
        var ret = '';

        if (isNaN(w) || isNaN(x) || isNaN(y) || isNaN(z)) {
            return 'NaN';
        }

        // Alternative design?
        // '(%f, [%f %f %f])'

        ret = numToStr(w, '', ret);
        ret += numToStr(x, 'i', ret);
        ret += numToStr(y, 'j', ret);
        ret += numToStr(z, 'k', ret);

        if ('' === ret)
            return '0';

        return ret;
    };

    /**
     * Returns the real part of the quaternion
     * @static
     * @param {Quaternion} quat
     * @return {*} 
     * @memberof Quaternion
     */
    static real(quat) {
        return quat.w;
    };

    /**
     * Returns the imaginary part of the quaternion as a 3D vector / array
     * @static
     * @param {Quaternion} quat
     * @return {*} 
     * @memberof Quaternion
     */
    static imag(quat) {
        return [0, quat.x, quat.y, quat.z];
    };

    /**
     * Gets the actual quaternion as a 4D vector / array
     * @static
     * @param {Quaternion} quat
     * @return {*} 
     * @memberof Quaternion
     */
    static toVector(/** @type {Quaternion} */ quat) {
        return [quat.w, quat.x, quat.y, quat.z];
    };

    /**
     * Calculates the 3x3 rotation matrix for the current quat
     * @static
     * @param {Quaternion} quat
     * @param {boolean} twoD
     * @return {*} 
     * @memberof Quaternion
     */
    static toMatrix(quat, twoD) {

        var w = quat.w;
        var x = quat.x;
        var y = quat.y;
        var z = quat.z;

        var wx = w * x, wy = w * y, wz = w * z;
        var xx = x * x, xy = x * y, xz = x * z;
        var yy = y * y, yz = y * z, zz = z * z;

        if (twoD) {
            return [
                [1 - 2 * (yy + zz), 2 * (xy - wz), 2 * (xz + wy)],
                [2 * (xy + wz), 1 - 2 * (xx + zz), 2 * (yz - wx)],
                [2 * (xz - wy), 2 * (yz + wx), 1 - 2 * (xx + yy)]];
        }

        return [
            1 - 2 * (yy + zz), 2 * (xy - wz), 2 * (xz + wy),
            2 * (xy + wz), 1 - 2 * (xx + zz), 2 * (yz - wx),
            2 * (xz - wy), 2 * (yz + wx), 1 - 2 * (xx + yy)
        ];
    };

    /**
     * Calculates the homogeneous 4x4 rotation matrix for the current quat
     * @static
     * @param {Quaternion} quat
     * @param {boolean} twoD
     * @return {*} 
     * @memberof Quaternion
     */
    static toMatrix4(quat, twoD) {

        var w = quat.w;
        var x = quat.x;
        var y = quat.y;
        var z = quat.z;

        var wx = w * x, wy = w * y, wz = w * z;
        var xx = x * x, xy = x * y, xz = x * z;
        var yy = y * y, yz = y * z, zz = z * z;

        if (twoD) {
            return [
                [1 - 2 * (yy + zz), 2 * (xy - wz), 2 * (xz + wy), 0],
                [2 * (xy + wz), 1 - 2 * (xx + zz), 2 * (yz - wx), 0],
                [2 * (xz - wy), 2 * (yz + wx), 1 - 2 * (xx + yy), 0],
                [0, 0, 0, 1]];
        }

        return [
            1 - 2 * (yy + zz), 2 * (xy - wz), 2 * (xz + wy), 0,
            2 * (xy + wz), 1 - 2 * (xx + zz), 2 * (yz - wx), 0,
            2 * (xz - wy), 2 * (yz + wx), 1 - 2 * (xx + yy), 0,
            0, 0, 0, 1];
    };

    /**
     * Calculates the Euler angles represented by the current quat
     * @static
     * @param {Quaternion} quat
     * @return {*} 
     * @memberof Quaternion
     */
    static toEuler(quat) {

        var w = quat.w;
        var x = quat.x;
        var y = quat.y;
        var z = quat.z;

        var t = 2 * (w * y - z * x);

        return {
            // X-axis rotation
            roll: Math.atan2(2 * (w * x + y * z), 1 - 2 * (x * x + y * y)),
            // Y-axis rotation
            pitch: t >= 1 ? Math.PI / 2 : (t <= -1 ? -Math.PI / 2 : Math.asin(t)),
            // Z-axis rotation
            yaw: Math.atan2(2 * (w * z + x * y), 1 - 2 * (y * y + z * z))
        };
    };

    /**
     * Clones the actual object
     * @static
     * @param {Quaternion} quat
     * @return {*} 
     * @memberof Quaternion
     */
    static clone(quat) {
        return new Quaternion(quat.w, quat.x, quat.y, quat.z);
    };


    /**
     * Rotates a vector according to the current quaternion, assumes |q|=1
     * @link https://www.xarg.org/proof/vector-rotation-using-quaternions/
     * @static
     * @param {Quaternion} quat
     * @param {*} v
     * @return {*} 
     * @memberof Quaternion
     */
    static rotateVector(quat, v) {

        var qw = quat.w;
        var qx = quat.x;
        var qy = quat.y;
        var qz = quat.z;

        var vx = v[0];
        var vy = v[1];
        var vz = v[2];

        // t = 2q x v
        var tx = 2 * (qy * vz - qz * vy);
        var ty = 2 * (qz * vx - qx * vz);
        var tz = 2 * (qx * vy - qy * vx);

        // v + w t + q x t
        return [
            vx + qw * tx + qy * tz - qz * ty,
            vy + qw * ty + qz * tx - qx * tz,
            vz + qw * tz + qx * ty - qy * tx
        ];
    };

    /**
     * Gets a function to spherically interpolate between two quaternions
     * @static
     * @param {Quaternion} quat1
     * @param {Quaternion} quat2
     * @return {*} 
     * @memberof Quaternion
     */
    static slerp(/** @type {Quaternion} */ quat1, /** @type {Quaternion} */ quat2) {
        var w1 = quat1.w;
        var x1 = quat1.x;
        var y1 = quat1.y;
        var z1 = quat1.z;

        var w2 = quat2.w;
        var x2 = quat2.x;
        var y2 = quat2.y;
        var z2 = quat2.z;

        var cosTheta0 = w1 * w2 + x1 * x2 + y1 * y2 + z1 * z2;

        if (cosTheta0 < 0) {
            w1 = -w1;
            x1 = -x1;
            y1 = -y1;
            z1 = -z1;
            cosTheta0 = -cosTheta0;
        }

        if (cosTheta0 >= 1 - Quaternion.epsilon) {
            return function (/** @type {*} */ pct) {
                return new Quaternion(
                    w1 + pct * (w2 - w1),
                    x1 + pct * (x2 - x1),
                    y1 + pct * (y2 - y1),
                    z1 + pct * (z2 - z1))['normalize']();
            };
        }

        var Theta0 = Math.acos(cosTheta0);
        var sinTheta0 = Math.sin(Theta0);

        return function (/** @type {*}*/ pct) {

            var Theta = Theta0 * pct;
            var sinTheta = Math.sin(Theta);
            var cosTheta = Math.cos(Theta);

            var s0 = cosTheta - cosTheta0 * sinTheta / sinTheta0;
            var s1 = sinTheta / sinTheta0;

            return new Quaternion(
                s0 * w1 + s1 * w2,
                s0 * x1 + s1 * x2,
                s0 * y1 + s1 * y2,
                s0 * z1 + s1 * z2);
        };
    };

    /**
     * @param {string} str
     */
    static fromString(str) {
        /**@type {Quaternion} */
        let resultQuat;
        str = str.replace(/\s/g, '');
        
        if (str.includes('i') || str.includes('j') || str.includes('k')) {
            const quaternion = { w: 0, x: 0, y: 0, z: 0 };

            const regex = /([+-]?\d*[ijk])|([+-]?\d+)/g;
            const matches = [...str.matchAll(regex)];

            for (const match of matches) {
                let value, type;

                // If the match includes an 'i', 'j', or 'k'
                if (match[1]) {
                    value = match[1].slice(0, -1);
                    if (value === '+' || value === '' || value === '-') {
                        value = value === '-' ? '-1' : '1';
                    }
                    type = match[1].slice(-1);
                } else {
                    value = match[2];
                    type = 'w';
                }
                
                value = parseInt(value);

                switch (type) {
                case 'i':
                    quaternion.x += value;
                    break;
                case 'j':
                    quaternion.y += value;
                    break;
                case 'k':
                    quaternion.z += value;
                    break;
                case 'w':
                    quaternion.w += value;
                    break;
                }
            }

            resultQuat = new Quaternion(quaternion.w, quaternion.x, quaternion.y, quaternion.z);

        } else {
            var digits = str.match(/\d+/g);
            /**@type {number[]} */
            var numbers = [];

            if (digits !== null) {
                numbers = digits.map(function (digit) {
                    return Number(digit);
                });
            }

            while (numbers.length < 4) {
                numbers.push(0);
            }

            resultQuat = new Quaternion(numbers[0], numbers[1], numbers[2], numbers[3]);
        }
        
        return resultQuat;
    }

    /**
     * @param {Quaternion} quat 
     */
    toQuatSignal(quat) {
        const Reactive = require('Reactive');
        return Reactive.quaternion(quat.w, quat.x, quat.y, quat.z);
    }
    
    /**
     * @param {QuaternionSignal} quat
     */
    fromQuatSignal(quat) {
        return new Quaternion(
            quat.w.pinLastValue(),
            quat.x.pinLastValue(),
            quat.y.pinLastValue(),
            quat.z.pinLastValue()
        )
    }

    add(/** @type {Quaternion} */ quat) {
        return Quaternion.add(this, quat);
    };

    sub(/** @type {Quaternion} */ quat) {
        return Quaternion.sub(this, quat);
    };

    neg() {
        return Quaternion.neg(this);
    };

    norm() {
        return Quaternion.norm(this);
    };

    normSq() {
        return Quaternion.normSq(this);
    };

    normalize() {
        return Quaternion.normalize(this);
    };

    mul(/** @type {Quaternion} */ quat) {
        return Quaternion.mul(this, quat);
    };

    scale(/** @type {number} */ s) {
        return Quaternion.scale(this, s);
    };

    dot(/** @type {Quaternion} */ quat) {
        return Quaternion.dot(this, quat);
    };

    inverse() {
        return Quaternion.inverse(this);
    };

    div(/** @type {Quaternion} */ quat) {
        return Quaternion.div(this, quat);
    };

    conjugate() {
        return Quaternion.conjugate(this);
    };

    exp() {
        return Quaternion.exp(this);
    };

    log() {
        return Quaternion.log(this);
    };

    pow(/** @type {Quaternion} */ quat) {
        return Quaternion.pow(this, quat);
    };

    equals(/** @type {Quaternion} */ quat) {
        return Quaternion.equals(this, quat);
    };

    isFinite() {
        return Quaternion.isFinite(this);
    };

    isNaN() {
        return Quaternion.isNaN(this);
    };

    toString() {
        return Quaternion.toString(this);
    };

    real() {
        return Quaternion.real(this);
    };

    imag() {
        return Quaternion.imag(this);
    };

    toVector() {
        return Quaternion.toVector(this);
    };

    toMatrix(/** @type {boolean} */ twoD) {
        return Quaternion.toMatrix(this, twoD);
    }

    toMatrix4(/** @type {boolean} */ twoD) {
        return Quaternion.toMatrix4(this, twoD);
    }

    toEuler() {
        return Quaternion.toEuler(this);
    };

    clone() {
        return Quaternion.clone(this);
    };

    rotateVector(/** @type {[]} */ v) {
        return Quaternion.rotateVector(this, v);
    };

    slerp(/** @type {Quaternion} */ quat) {
        return Quaternion.slerp(this, quat);
    };
}
/**
 * Calculates log(sqrt(a^2+b^2)) in a way to avoid overflows
 *
 * @param {number} a
 * @param {number} b
 * @returns {number}
 */
function logHypot(a, b) {

    var _a = Math.abs(a);
    var _b = Math.abs(b);

    if (a === 0) {
        return Math.log(_b);
    }

    if (b === 0) {
        return Math.log(_a);
    }

    if (_a < 3000 && _b < 3000) {
        return 0.5 * Math.log(a * a + b * b);
    }

    a = a / 2;
    b = b / 2;

    return 0.5 * Math.log(a * a + b * b) + Math.LN2;
}

/**
 *
 * @param {*} n
 * @param {*} char
 * @param {*} prev
 * @return {*} 
 */
function numToStr(n, char, prev) {

    var ret = '';

    if (n !== 0) {

        if (prev !== '') {
            ret += n < 0 ? ' - ' : ' + ';
        } else if (n < 0) {
            ret += '-';
        }

        n = Math.abs(n);

        if (1 !== n || char === '') {
            ret += n;
        }
        ret += char;
    }
    return ret;
}

/**
 * Class containing auxiliary mathematical functions
 *
 * @export
 * @class AMath
 */
class AMath {
	/**
	 * Converts degrees to radians
	 *
	 * @export
	 * @param {number} degrees
	 * @returns {number}
	 */
	static Deg2Rad(degrees) {
		return degrees * (Math.PI / 180);
	}

	/**
	 * Converts radians to degrees
	 *
	 * @export
	 * @param {number} radians
	 * @returns {number}
	 */
	static Rad2Deg(radians) {
		return radians * (180 / Math.PI);
	}

	/**
	 * Returns a random integer in the range of numbers [min, max)
	 *
	 * @export
	 * @param {number} lowerBound
	 * @param {number} upperBound
	 * @returns {number}
	 */
	static GetRandomInt(lowerBound, upperBound) {
		return Math.floor(Math.random() * ((upperBound) - lowerBound) + lowerBound);
	}

	/**
	 * Returns a random float in the range of numbers [min, max)
	 *
	 * @export
	 * @param {number} lowerBound
	 * @param {number} upperBound
	 * @returns {number}
	 */
	static GetRandomFloat(lowerBound, upperBound) {
		return Math.random() * ((upperBound) - lowerBound) + lowerBound;
	}

	/**
	 * Checks if the given number is in the range of numbers [lowerBound, upperBound]
	 *
	 * @export
	 * @param {number} number
	 * @param {number} lowerBound
	 * @param {number} upperBound
	 * @returns {boolean}
	 */
	static InRange(number, lowerBound, upperBound) {
		return ((number - lowerBound) * (number - upperBound) <= 0);
	}

	/**
	 * Returns a random item based on weights or item drop chances
	 *
	 * @static
	 * @param {Object.<string, number>} weightedArray Object where keys are elements and values are weights as a number
	 * @returns {string}
	 * @memberof AMath
	 */
	static chooseWeighted(weightedArray) {
		let items = Object.keys(weightedArray);
		let chances = Object.values(weightedArray);

		// Calculate the sum of all weights
		let sum = chances.reduce((total, element) => total + element, 0);
		let total = 0;
		chances = chances.map(chance => total += chance);
		let rand = Math.random() * sum;
		let result = items[chances.filter(el => el <= rand).length];
		if (result === '0') return null;
		else return result;
	}


	/**
	 * Returns a random element of the given array
	 *
	 * @static
	 * @param {any[]} array
	 * @returns {any}
	 * @memberof AMath
	 */
	static chooseRandom(array) {
		return array[this.GetRandomInt(0, array.length)];
	}

	/**
	 * Get distance between coordinates
	 *
	 * @static
	 * @param {number} x
	 * @param {number} y
	 * @returns {number}
	 * @memberof AMath
	 */
	static GetDistance(x, y) {
		return Math.sqrt((x - y) * (x - y))
	}

	/**
	 * Returns the interpolated value between start and end
	 *
	 * @static
	 * @param {number} start Start value
	 * @param {number} end End value
	 * @param {number} t Interpolation (from 0 to 1)
	 * @returns {number}
	 * @memberof AMath
	 */
	static lerp(start, end, t) {
		return start * (1 - t) + end * t;
	}

	/**
	 * Checks if two line segments intersect
	 *
	 * @static
	 * @param {number} positionA
	 * @param {number} lengthA
	 * @param {number} positionB
	 * @param {number} lengthB
	 * @returns {boolean}
	 * @memberof AMath
	 */
	static checkLinesIntersection(positionA, lengthA, positionB, lengthB) {
		return Math.abs(positionA - positionB) <= ((lengthA / 2) + (lengthB / 2));
	}

	/**
	 *
	 *
	 * @static
	 * @param {Vector3} point
	 * @param {Vector3} rectCenter
	 * @param {Vector3} rectSize
	 * @param {Quaternion} rectRotation
	 * @memberof AMath
	 */
	static isPointInRect(point, rectCenter, rectSize, rectRotation) {
		let pointInRectFrame = Vector3.rotate(
			point.sub(rectCenter),
			rectRotation.inverse()
		);

		return Math.abs(pointInRectFrame.x) <= (rectSize.x / 2) && Math.abs(pointInRectFrame.z) <= (rectSize.y / 2);
	}

	/**
	 *
	 *
	 * @static
	 * @param {Vector3} point
	 * @param {Vector3} boxCenter
	 * @param {Vector3} boxSize
	 * @param {Quaternion} boxRotation
	 * @memberof AMath
	 */
	static isPointInBox(point, boxCenter, boxSize, boxRotation) {
		let pointInBoxFrame = Vector3.rotate(
			point.sub(boxCenter),
			boxRotation.inverse()
		);

		return (
			Math.abs(pointInBoxFrame.x) <= (boxSize.x / 2) &&
			Math.abs(pointInBoxFrame.y) <= (boxSize.y / 2) &&
			Math.abs(pointInBoxFrame.z) <= (boxSize.z / 2)
		);
	}

	/**
	 * @static
	 * @param {number} value
	 * @param {number} min
	 * @param {number} max
	 * @returns {number}
	 * @memberof AMath
	 */
	static clamp(value, min, max) {
		return Math.min(Math.max(value, min), max);
	}
}

var AMath$1 = AMath;

class Matrix {
    constructor(elements = [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]) {
        this.elements = elements;
    }

    /**
     * @param {Matrix} matrix
     * @memberof Matrix
     */
    multiply(matrix) {
        // Multiply this matrix by the given matrix
        const elements = [];
        for (let i = 0; i < 16; i++) {
            elements[i] = 0;
            for (let j = 0; j < 4; j++) {
                elements[i] += this.elements[j * 4 + i % 4] * matrix.elements[j + Math.floor(i / 4) * 4];
            }
        }

        return new Matrix(elements)
    }

    /**
     * @param {Vector3} point
     * @return {Vector3} 
     * @memberof Matrix
     */
    transformPoint(point) {
        const elements = this.elements;

        const transformedPoint = Vector3.fromArray([
            point.x * elements[0] + point.y * elements[4] + point.z * elements[8] + elements[12],
            point.x * elements[1] + point.y * elements[5] + point.z * elements[9] + elements[13],
            point.x * elements[2] + point.y * elements[6] + point.z * elements[10] + elements[14],
        ]);

        return transformedPoint;
    }

    /**
     * @param {Quaternion} quaternion
     * @return {Quaternion} 
     * @memberof Matrix
     */
    transformQuaternion(quaternion) {
        const x = quaternion.x;
        const y = quaternion.y;
        const z = quaternion.z;
        const w = quaternion.w;

        const elements = this.elements;

        const array = [
            x * elements[0] + y * elements[4] + z * elements[8] + w * elements[12],
            x * elements[1] + y * elements[5] + z * elements[9] + w * elements[13],
            x * elements[2] + y * elements[6] + z * elements[10] + w * elements[14],
            x * elements[3] + y * elements[7] + z * elements[11] + w * elements[15]
        ];

        const transformedQuaternion = new Quaternion(
            array[3],
            array[0],
            array[1],
            array[2],
        );

        return transformedQuaternion;
    }

    /**
     * @param {Vector3} vector
     * @return {Vector3} 
     * @memberof Matrix
     */
    transformVector(vector) {
        const x = vector.x;
        const y = vector.y;
        const z = vector.z;

        const elements = this.elements;

        const transformedVector = Vector3.fromArray([
            x * elements[0] + y * elements[4] + z * elements[8],
            x * elements[1] + y * elements[5] + z * elements[9],
            x * elements[2] + y * elements[6] + z * elements[10]
        ]);

        return transformedVector;
    }

    inverse() {
        const te = this.elements;
        const n11 = te[0], n21 = te[1], n31 = te[2], n41 = te[3];
        const n12 = te[4], n22 = te[5], n32 = te[6], n42 = te[7];
        const n13 = te[8], n23 = te[9], n33 = te[10], n43 = te[11];
        const n14 = te[12], n24 = te[13], n34 = te[14], n44 = te[15];

        const t11 = n23 * n34 * n42 - n24 * n33 * n42 + n24 * n32 * n43 - n22 * n34 * n43 - n23 * n32 * n44 + n22 * n33 * n44;
        const t12 = n14 * n33 * n42 - n13 * n34 * n42 - n14 * n32 * n43 + n12 * n34 * n43 + n13 * n32 * n44 - n12 * n33 * n44;
        const t13 = n13 * n24 * n42 - n14 * n23 * n42 + n14 * n22 * n43 - n12 * n24 * n43 - n13 * n22 * n44 + n12 * n23 * n44;
        const t14 = n14 * n23 * n32 - n13 * n24 * n32 - n14 * n22 * n33 + n12 * n24 * n33 + n13 * n22 * n34 - n12 * n23 * n34;

        const det = n11 * t11 + n21 * t12 + n31 * t13 + n41 * t14;

        if (det === 0) throw new Error('Det === 0')

        const detInv = 1 / det;

        te[0] = t11 * detInv;
        te[1] = (n24 * n33 * n41 - n23 * n34 * n41 - n24 * n31 * n43 + n21 * n34 * n43 + n23 * n31 * n44 - n21 * n33 * n44) * detInv;
        te[2] = (n22 * n34 * n41 - n24 * n32 * n41 + n24 * n31 * n42 - n21 * n34 * n42 - n22 * n31 * n44 + n21 * n32 * n44) * detInv;
        te[3] = (n23 * n32 * n41 - n22 * n33 * n41 - n23 * n31 * n42 + n21 * n33 * n42 + n22 * n31 * n43 - n21 * n32 * n43) * detInv;

        te[4] = t12 * detInv;
        te[5] = (n13 * n34 * n41 - n14 * n33 * n41 + n14 * n31 * n43 - n11 * n34 * n43 - n13 * n31 * n44 + n11 * n33 * n44) * detInv;
        te[6] = (n14 * n32 * n41 - n12 * n34 * n41 - n14 * n31 * n42 + n11 * n34 * n42 + n12 * n31 * n44 - n11 * n32 * n44) * detInv;
        te[7] = (n12 * n33 * n41 - n13 * n32 * n41 + n13 * n31 * n42 - n11 * n33 * n42 - n12 * n31 * n43 + n11 * n32 * n43) * detInv;

        te[8] = t13 * detInv;
        te[9] = (n14 * n23 * n41 - n13 * n24 * n41 - n14 * n21 * n43 + n11 * n24 * n43 + n13 * n21 * n44 - n11 * n23 * n44) * detInv;
        te[10] = (n12 * n24 * n41 - n14 * n22 * n41 + n14 * n21 * n42 - n11 * n24 * n42 - n12 * n21 * n44 + n11 * n22 * n44) * detInv;
        te[11] = (n13 * n22 * n41 - n12 * n23 * n41 - n13 * n21 * n42 + n11 * n23 * n42 + n12 * n21 * n43 - n11 * n22 * n43) * detInv;

        te[12] = t14 * detInv;
        te[13] = (n13 * n24 * n31 - n14 * n23 * n31 + n14 * n21 * n33 - n11 * n24 * n33 - n13 * n21 * n34 + n11 * n23 * n34) * detInv;
        te[14] = (n14 * n22 * n31 - n12 * n24 * n31 - n14 * n21 * n32 + n11 * n24 * n32 + n12 * n21 * n34 - n11 * n22 * n34) * detInv;
        te[15] = (n12 * n23 * n31 - n13 * n22 * n31 + n13 * n21 * n32 - n11 * n23 * n32 - n12 * n21 * n33 + n11 * n22 * n33) * detInv;

        return new Matrix(te);

    }

    /**
     * @param {Vector3} target
     * @param {Vector3} [up=new Vector3(0, 1, 0)]
     * @return {Matrix} 
     * @memberof Matrix
     */
    lookAt(target, up = new Vector3(0, 1, 0)) {
        let z = new Vector3(this.elements[12], this.elements[13], this.elements[14]).sub(target).normalize();
        let x = up.cross(z).normalize();
        let y = z.cross(x);

        const te = this.elements;

        const scale = new Vector3(
            Math.sqrt(te[0] * te[0] + te[1] * te[1] + te[2] * te[2]),
            Math.sqrt(te[4] * te[4] + te[5] * te[5] + te[6] * te[6]),
            Math.sqrt(te[8] * te[8] + te[9] * te[9] + te[10] * te[10])
        );

        let el = [];
        el[0] = x.x * scale.x;
        el[1] = y.x * scale.x;
        el[2] = z.x * scale.x;
        el[3] = 0;
        el[4] = x.y * scale.y;
        el[5] = y.y * scale.y;
        el[6] = z.y * scale.y;
        el[7] = 0;
        el[8] = x.z * scale.z;
        el[9] = y.z * scale.z;
        el[10] = z.z * scale.z;
        el[11] = 0;
        el[12] = this.elements[12];
        el[13] = this.elements[13];
        el[14] = this.elements[14];
        el[15] = 1;

        return new Matrix(el);
    }

    transpose() {
        const transposed = [];
        for (let i = 0; i < 4; i++) {
            for (let j = 0; j < 4; j++) {
                transposed.push(this.elements[j * 4 + i]);
            }
        }

        return new Matrix(transposed);
    }

    /**
     * @static
     * @param {Vector3} position
     * @param {Quaternion} quaternion
     * @param {Vector3} scale
     * @return {Matrix} 
     * @memberof Matrix
     */
    static compose(position, quaternion, scale) {
        const te = [];

        const x = quaternion.x, y = quaternion.y, z = quaternion.z, w = quaternion.w;
        const x2 = x + x, y2 = y + y, z2 = z + z;
        const xx = x * x2, xy = x * y2, xz = x * z2;
        const yy = y * y2, yz = y * z2, zz = z * z2;
        const wx = w * x2, wy = w * y2, wz = w * z2;

        const sx = scale.x, sy = scale.y, sz = scale.z;

        te[0] = (1 - (yy + zz)) * sx;
        te[1] = (xy + wz) * sx;
        te[2] = (xz - wy) * sx;
        te[3] = 0;

        te[4] = (xy - wz) * sy;
        te[5] = (1 - (xx + zz)) * sy;
        te[6] = (yz + wx) * sy;
        te[7] = 0;

        te[8] = (xz + wy) * sz;
        te[9] = (yz - wx) * sz;
        te[10] = (1 - (xx + yy)) * sz;
        te[11] = 0;

        te[12] = position.x;
        te[13] = position.y;
        te[14] = position.z;
        te[15] = 1;

        return new Matrix(te);
    }

    toString() {
        let str = '';

        for (let i = 0; i < 4; i++) {
            for (let j = 0; j < 4; j++) {
                str += this.elements[i * 4 + j] + ' ';
            }
    
            str += '\n';
        }

        return str
    }
}

class Vector2 {
	/**
	 * Creates an instance of Vector3.
	 * @param {number} x
	 * @param {number} y
	 * @memberof Vector2
	 */
	constructor(x, y) {
		Object.defineProperty(this, "x", { get: () => x });
		Object.defineProperty(this, "y", { get: () => y });
	}

	get x() { return 0; }
	get y() { return 0; }

	static get zero() { return new Vector2(0, 0); }
	static get right() { return new Vector2(1, 0); }
	static get up() { return new Vector2(0, 1); }

	get magnitude() {
		return Math.sqrt(this.x ** 2 + this.y ** 2);
	}

	get magnitudeSquared() {
		return this.x ** 2 + this.y ** 2;
	}

	static add(/** @type {Vector2} */ vec1, /** @type {Vector2} */ vec2) {
		return new Vector2(
			vec1.x + vec2.x,
			vec1.y + vec2.y
		);
	};

	static sub(/** @type {Vector2} */ vec1, /** @type {Vector2} */ vec2) {
		return new Vector2(
			vec1.x - vec2.x,
			vec1.y - vec2.y
		);
	};

	static mul(/** @type {Vector2} */ vec1, /** @type {number} */ scalar) {
		return new Vector2(
			vec1.x * scalar,
			vec1.y * scalar
		);
	};

	static div(/** @type {Vector2} */ vec1, /** @type {number} */ scalar) {
		return new Vector2(
			vec1.x / scalar,
			vec1.y / scalar
		);
	};

	/**
	 * Returns dot product of two vectors
	 *
	 * @static
	 * @param {Vector2} vec1
	 * @param {Vector2} vec2
	 * @returns {number}
	 * @memberof Vector2
	 */
	static dot(vec1, vec2) {
		return vec1.x * vec2.x + vec1.y * vec2.y;
	};

	/**
	 * Returns cross product of two vectors
	 *
	 * @static
	 * @param {Vector2} vec1
	 * @param {Vector2} vec2
	 * @returns {number}
	 * @memberof Vector2
	 */
	static cross(vec1, vec2) {
		return vec1.x * vec2.y - vec1.y * vec2.x;
	};

	/**
	 * Returns the angle between two vectors in radians
	 *
	 * @static
	 * @param {Vector2} vec1
	 * @param {Vector2} vec2
	 * @returns
	 * @memberof Vector2
	 */
	static angleBetween(vec1, vec2) {
		let cos = Vector2.dot(vec1, vec2) / (vec1.magnitude * vec2.magnitude);
		return Math.acos(cos > 1 || cos < -1 ? Math.round(cos) : cos);
	}

	/**
	 * Returns the number which is the distance between two vectors
	 * @param {Vector2} vector1
	 * @param {Vector2} vector2 
	 */
	static distance(vector1, vector2) {
		let subVec = new Vector2(
			vector1.x - vector2.x,
			vector1.y - vector2.y
		);

		return Math.sqrt(subVec.x ** 2 + subVec.y ** 2);
	}

	/**
	 * Linearly interpolates between two vectors
	 *
	 * @param {Vector2} vec1
	 * @param {Vector2} vec2
	 * @param {number} t Interpolation factor - should be in range [0; 1]
	 * @return {Vector2} 
	 */
	static lerp(vec1, vec2, t) {
		t = AMath$1.clamp(t, 0, 1);
		return new Vector2(
			vec1.x + (vec2.x - vec1.x) * t,
			vec1.y + (vec2.y - vec1.y) * t
		);
	}

	add(/** @type {Vector2} */ vec) {
		return Vector2.add(this, vec);
	};

	sub(/** @type {Vector2} */ vec) {
		return Vector2.sub(this, vec);
	};

	/**
	 *
	 *
	 * @param {number} scalar
	 * @memberof Vector2
	 */
	mul(scalar) {
		return Vector2.mul(this, scalar);
	}

	/**
	 *
	 *
	 * @param {number} scalar
	 * @memberof Vector2
	 */
	div(scalar) {
		return Vector2.div(this, scalar);
	}

	/**
	 * Returns dot product of two vectors
	 *
	 * @static
	 * @param {Vector2} vec
	 * @returns {number}
	 * @memberof Vector2
	 */
	dot(vec) {
		return Vector2.dot(this, vec);
	};

	/**
	 * Returns cross product of two vectors
	 *
	 * @param {Vector2} vec
	 * @returns {number}
	 * @memberof Vector2
	 */
	cross(vec) {
		return Vector2.cross(this, vec);
	};

	/**
	 * Returns the angle between two vectors in radians
	 *
	 * @static
	 * @param {Vector2} vec
	 * @returns
	 * @memberof Vector2
	 */
	angleBetween(vec) {
		return Vector2.angleBetween(this, vec);
	}

	/**
	 * Returns the number which is the distance between two vectors
	 * @param {Vector2} vector 
	 */
	distance(vector) {
		return Vector2.distance(this, vector);
	}

	/**
	 * Returns the vector, the madnitude of which is limited by specified one
	 *
	 * @param {number} maxLength
	 * @returns
	 */
	clampMagnitude(maxLength) {
		let sqrmag = this.magnitudeSquared;
		if (sqrmag > maxLength * maxLength) {
			let mag = Math.sqrt(sqrmag);
			let normalized_x = this.x / mag;
			let normalized_y = this.y / mag;
			return new Vector2(
				normalized_x * maxLength,
				normalized_y * maxLength,
			);
		}
		return this;
	}

	/**
	 * Linearly interpolates between two vectors
	 *
	 * @param {Vector2} vector
	 * @param {number} t Interpolation factor - should be in range [0; 1]
	 * @return {Vector2} 
	 */
	lerp(vector, t) {
		return Vector2.lerp(this, vector, t);
	}

	neg() {
		return this.mul(-1);
	}

	normalize() {
		return this.div(this.magnitude);
	}

	/**
	 * @static
	 * @param {number[]} array
	 * @returns
	 * @memberof Vector2
	 */
	static fromArray(array) {
		return new Vector2(
			array[0],
			array[1]
		)
	}

	toPoint2dSignal() {
		const Reactive = require('Reactive');
		Reactive.point2d(this.x, this.y);
	}

	toString(digits = 2) {
		return "(" + this.x.toFixed(digits) + "; " + this.y.toFixed(digits) + ")";
	}

	/**
	 * @return {[number, number]} 
	 * @memberof Vector2
	 */
	toArray() {
		return [this.x, this.y];
	}
}

export { AMath, Matrix, Quaternion, Vector2, Vector3 };
