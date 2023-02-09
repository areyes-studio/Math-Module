import Vector3 from "./Vector3";
import Quaternion from "quaternion";

/**
 * Class containing auxiliary mathematical functions
 *
 * @export
 * @class AMath
 */
export class AMath {
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

export default AMath;
