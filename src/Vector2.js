import AMath from "./AMath";

export default class Vector2 {
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
		let cos = Vector2.dot(vec1, vec2) / (vec1.magnitude * vec2.magnitude)
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
		)

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
		t = AMath.clamp(t, 0, 1);
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