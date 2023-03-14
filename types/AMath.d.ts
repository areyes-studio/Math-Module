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
    static Deg2Rad(degrees: number): number;
    /**
     * Converts radians to degrees
     *
     * @export
     * @param {number} radians
     * @returns {number}
     */
    static Rad2Deg(radians: number): number;
    /**
     * Returns a random integer in the range of numbers [min, max)
     *
     * @export
     * @param {number} lowerBound
     * @param {number} upperBound
     * @returns {number}
     */
    static GetRandomInt(lowerBound: number, upperBound: number): number;
    /**
     * Returns a random float in the range of numbers [min, max)
     *
     * @export
     * @param {number} lowerBound
     * @param {number} upperBound
     * @returns {number}
     */
    static GetRandomFloat(lowerBound: number, upperBound: number): number;
    /**
     * Checks if the given number is in the range of numbers [lowerBound, upperBound]
     *
     * @export
     * @param {number} number
     * @param {number} lowerBound
     * @param {number} upperBound
     * @returns {boolean}
     */
    static InRange(number: number, lowerBound: number, upperBound: number): boolean;
    /**
     * Returns a random item based on weights or item drop chances
     *
     * @static
     * @param {Object.<string, number>} weightedArray Object where keys are elements and values are weights as a number
     * @returns {string}
     * @memberof AMath
     */
    static chooseWeighted(weightedArray: {
        [x: string]: number;
    }): string;
    /**
     * Returns a random element of the given array
     *
     * @static
     * @param {any[]} array
     * @returns {any}
     * @memberof AMath
     */
    static chooseRandom(array: any[]): any;
    /**
     * Get distance between coordinates
     *
     * @static
     * @param {number} x
     * @param {number} y
     * @returns {number}
     * @memberof AMath
     */
    static GetDistance(x: number, y: number): number;
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
    static lerp(start: number, end: number, t: number): number;
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
    static checkLinesIntersection(positionA: number, lengthA: number, positionB: number, lengthB: number): boolean;
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
    static isPointInRect(point: Vector3, rectCenter: Vector3, rectSize: Vector3, rectRotation: Quaternion): boolean;
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
    static isPointInBox(point: Vector3, boxCenter: Vector3, boxSize: Vector3, boxRotation: Quaternion): boolean;
    /**
     * @static
     * @param {number} value
     * @param {number} min
     * @param {number} max
     * @returns {number}
     * @memberof AMath
     */
    static clamp(value: number, min: number, max: number): number;
}
export default AMath;
import Vector3 from "./Vector3.js";
import Quaternion from "quaternion";
