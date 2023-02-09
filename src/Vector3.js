import AMath from "./AMath";

export default class Vector3 {
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
		)

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
		t = AMath.clamp(t, 0, 1);
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
