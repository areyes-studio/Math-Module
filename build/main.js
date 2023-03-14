'use strict';

var Quaternion = require('quaternion');

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
	 * Возвращает скалярное произведение двух векторов (dot product)
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
	 * Возвращает векторное произведение двух векторов (cross product)
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
	 * Возвращает угол между двумя векторами в радианах
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
	 * Вращает вектор с помощью кватерниона
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
	 * Возвращает скалярное произведение двух векторов (dot product)
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
	 * Возвращает векторное произведение двух векторов (cross product)
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
	 * Возвращает угол между двумя векторами в радианах
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
	 * Возвращает вектор, длина которого ограничена заданной
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
	 * Возвращает скалярное произведение двух векторов (dot product)
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
	 * Возвращает вектор, длина которого ограничена заданной
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

Object.defineProperty(exports, 'Quaternion', {
	enumerable: true,
	get: function () { return Quaternion.Quaternion; }
});
exports.AMath = AMath;
exports.Matrix = Matrix;
exports.Vector2 = Vector2;
exports.Vector3 = Vector3;
