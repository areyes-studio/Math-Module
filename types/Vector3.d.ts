export class Vector3 {
    static get zero(): Vector3;
    static get right(): Vector3;
    static get up(): Vector3;
    static get forward(): Vector3;
    static add(vec1: Vector3, vec2: Vector3): Vector3;
    static sub(vec1: Vector3, vec2: Vector3): Vector3;
    static mul(vec1: Vector3, scalar: number): Vector3;
    static div(vec1: Vector3, scalar: number): Vector3;
    /**
     * Returns dot product of two vectors
     *
     * @static
     * @param {Vector3} vec1
     * @param {Vector3} vec2
     * @returns {number}
     * @memberof Vector3
     */
    static dot(vec1: Vector3, vec2: Vector3): number;
    /**
     * Returns cross product of two vectors
     *
     * @static
     * @param {Vector3} vec1
     * @param {Vector3} vec2
     * @returns
     * @memberof Vector3
     */
    static cross(vec1: Vector3, vec2: Vector3): Vector3;
    /**
     * Returns an angle between two vectors in radians
     *
     * @static
     * @param {Vector3} vec1
     * @param {Vector3} vec2
     * @returns
     * @memberof Vector3
     */
    static angleBetween(vec1: Vector3, vec2: Vector3): number;
    /**
     * Returns signed angle between two vectors on specified axis
     *
     * @static
     * @param {Vector3} vec1
     * @param {Vector3} vec2
     * @param {Vector3} axis
     * @memberof Vector3
     */
    static angleBetweenSigned(vec1: Vector3, vec2: Vector3, axis: Vector3): number;
    /**
     * Returns the number which is the distance between two vectors
     * @param {Vector3} vector1
     * @param {Vector3} vector2
     */
    static distance(vector1: Vector3, vector2: Vector3): number;
    /**
     * Rotates vector by quaternion
     *
     * @static
     * @param {*} v
     * @param {*} q
     * @returns
     * @memberof Vector3
     */
    static rotate(v: Vector3, q: {
        w: number;
        x: number;
        y: number;
        z: number;
    }): Vector3;
    /**
     * Projects a vector onto a plane defined by normal
     *
     * @static
     * @param {Vector3} vector
     * @param {Vector3} planeNormal
     * @memberof Vector3
     */
    static projectOnPlane(vector: Vector3, planeNormal: Vector3): Vector3;
    /**
     * Linearly interpolates between two vectors
     *
     * @param {Vector3} vec1
     * @param {Vector3} vec2
     * @param {number} t Interpolation factor - should be in range [0; 1]
     * @return {Vector3}
     */
    static lerp(vec1: Vector3, vec2: Vector3, t: number): Vector3;
    /**
     * @static
     * @param {PointSignal} point
     * @returns
     * @memberof Vector3
     */
    static fromPointSignal(point: PointSignal): Vector3;
    /**
     *
     *
     * @static
     * @param {number[]} array
     * @returns
     * @memberof Vector3
     */
    static fromArray(array: number[]): Vector3;
    /**
     * Creates an instance of Vector3.
     * @param {number} x
     * @param {number} y
     * @param {number} z
     * @memberof Vector3
     */
    constructor(x: number, y: number, z: number);
    get x(): number;
    get y(): number;
    get z(): number;
    add(vec: Vector3): Vector3;
    sub(vec: Vector3): Vector3;
    /**
     *
     *
     * @param {number} scalar
     * @memberof Vector3
     */
    mul(scalar: number): Vector3;
    /**
     *
     *
     * @param {number} scalar
     * @memberof Vector3
     */
    div(scalar: number): Vector3;
    /**
     * Returns dot product of two vectors
     *
     * @static
     * @param {Vector3} vec
     * @returns {number}
     * @memberof Vector3
     */
    dot(vec: Vector3): number;
    /**
     * Returns cross product of two vectors
     *
     * @static
     * @param {Vector3} vec
     * @returns
     * @memberof Vector3
     */
    cross(vec: Vector3): Vector3;
    /**
     * Returns an angle between two vectors in radians
     *
     * @static
     * @param {Vector3} vec
     * @returns
     * @memberof Vector3
     */
    angleBetween(vec: Vector3): number;
    /**
     * Returns signed angle between two vectors on specified axis
     *
     * @static
     * @param {Vector3} vector
     * @param {Vector3} axis
     * @memberof Vector3
     */
    angleBetweenSigned(vector: Vector3, axis: Vector3): number;
    /**
     * Returns the number which is the distance between two vectors
     * @param {Vector3} vector
     */
    distance(vector: Vector3): number;
    rotate(q: {
        w: number;
        x: number;
        y: number;
        z: number;
    }): Vector3;
    /**
     * Returns the vector, the madnitude of which is limited by specified one
     *
     * @param {number} maxLength
     * @returns
     */
    clampMagnitude(maxLength: number): Vector3;
    /**
     * Projects a vector onto a plane defined by normal
     *
     * @param {Vector3} planeNormal
     * @memberof Vector3
     */
    projectOnPlane(planeNormal: Vector3): Vector3;
    /**
     * Linearly interpolates between two vectors
     *
     * @param {Vector3} vector
     * @param {number} t Interpolation factor - should be in range [0; 1]
     * @return {Vector3}
     */
    lerp(vector: Vector3, t: number): Vector3;
    neg(): Vector3;
    get magnitude(): number;
    get magnitudeSquared(): number;
    normalize(): Vector3;
    toPointSignal(): PointSignal;
    toVectorSignal(): VectorSignal;
    toString(digits?: number): string;
    /**
     * @return {[number, number, number]}
     * @memberof Vector3
     */
    toArray(): [number, number, number];
}
export default Vector3;
