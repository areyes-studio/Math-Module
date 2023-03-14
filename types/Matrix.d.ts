export class Matrix {
    /**
     * @static
     * @param {Vector3} position
     * @param {Quaternion} quaternion
     * @param {Vector3} scale
     * @return {Matrix}
     * @memberof Matrix
     */
    static compose(position: Vector3, quaternion: Quaternion, scale: Vector3): Matrix;
    constructor(elements?: number[]);
    elements: number[];
    /**
     * @param {Matrix} matrix
     * @memberof Matrix
     */
    multiply(matrix: Matrix): Matrix;
    /**
     * @param {Vector3} point
     * @return {Vector3}
     * @memberof Matrix
     */
    transformPoint(point: Vector3): Vector3;
    /**
     * @param {Quaternion} quaternion
     * @return {Quaternion}
     * @memberof Matrix
     */
    transformQuaternion(quaternion: Quaternion): Quaternion;
    /**
     * @param {Vector3} vector
     * @return {Vector3}
     * @memberof Matrix
     */
    transformVector(vector: Vector3): Vector3;
    inverse(): Matrix;
    /**
     * @param {Vector3} target
     * @param {Vector3} [up=new Vector3(0, 1, 0)]
     * @return {Matrix}
     * @memberof Matrix
     */
    lookAt(target: Vector3, up?: Vector3): Matrix;
    transpose(): Matrix;
    toString(): string;
}
export default Matrix;
import Vector3 from "./Vector3.js";
import Quaternion from "quaternion";
