
import Quaternion from "quaternion";
import AMath from "../src/AMath";
import Vector3 from "../src/Vector3";

test('Converts degrees to radians - Deg2Rad method', () => {
    expect(AMath.Deg2Rad(180)).toBe(Math.PI);
    expect(AMath.Deg2Rad(0)).toBe(0);
    expect(AMath.Deg2Rad(45)).toBe(Math.PI / 4);
})

test('Convert radians to degrees - Rad2Deg method', () => {
    expect(AMath.Rad2Deg(Math.PI)).toBe(180);
    expect(AMath.Rad2Deg(0)).toBe(0);
    expect(AMath.Rad2Deg(Math.PI / 4)).toBe(45);
})

test('Returns a random integer in the range of numbers [min, max) - GetRandomInt method', () => {
    expect(AMath.GetRandomInt(-9, 10)).toBeGreaterThanOrEqual(-9);
    expect(AMath.GetRandomInt(-9, 10)).toBeLessThan(10)
    expect(AMath.GetRandomInt(10, 10)).toBe(10)
})

test('Returns a random number in the range of numbers [min, max)- GetRandomFloat method', () => {
    expect(AMath.GetRandomFloat(-9, 10)).toBeGreaterThanOrEqual(-9);
    expect(AMath.GetRandomFloat(-9, 10)).toBeLessThan(10)
    expect(AMath.GetRandomFloat(10, 10)).toBe(10)
})

test('Checks if the given number is in the range of numbers [lowerBound, upperBound] - InRange method', () => {
    expect(AMath.InRange(6, 6, 10)).toBe(true);
    expect(AMath.InRange(10, 6, 10)).toBe(true);
    expect(AMath.InRange(-7, 6, 10)).toBe(false);
})

test('Returns a random item based on the weight or drop chance of the item - chooseWeighted method', () => {
    let obj = {
        'a': 0,
        'b': 10,
        'c': 0
    }
    expect(AMath.chooseWeighted(obj)).toBe('b')
})

test("Get distance between coordinates (dimension 1) - GetDistance method", () => {
    expect(AMath.GetDistance(6.99999, 0.99999)).toBe(6)
    expect(AMath.GetDistance(-10, 20)).toBe(30)
    expect(AMath.GetDistance(-0, 0)).toBe(0)
})

test('Returns the interpolated value between the start and end - lerp method', () => {
    expect(AMath.lerp(10, 20, 0.5)).toBe(15);
    expect(AMath.lerp(1, 200, 0)).toBe(1);
    expect(AMath.lerp(-1.125, 200, 1)).toBe(200);

})

test('Checks if two line segments intersect on a line one-dimensional case - checkLinesIntersection method', () => {
    expect(AMath.checkLinesIntersection(0, 5, 10, 4)).toBe(false)
    expect(AMath.checkLinesIntersection(0, 5, 10, 16)).toBe(true)
    expect(AMath.checkLinesIntersection(-10, 40, 10, 0)).toBe(true)
})

test('reurn value or min or max when out of range[min, max] - clamp mehod', () => {
    expect(AMath.clamp(-9, 2, 3)).toBe(2)
    expect(AMath.clamp(10, 10, 10)).toBe(10)
    expect(AMath.clamp(-526.235, -526.1, -526.2)).toBe(-526.2)
})

test('determines whether a point belongs to a parallelepiped or not - isPointInBOX method', () => {
    expect(AMath.isPointInBox(new Vector3(0.5, 1, 2), new Vector3(0, 0, 0), new Vector3(1, 2, 4), new Quaternion(1, 1, 0, 0))).toBe(true)
    expect(AMath.isPointInBox(new Vector3(0, 0, 0), new Vector3(0, 0, 0), new Vector3(0, 0, 0), new Quaternion(1, 1, 0, 0))).toBe(true)
    expect(AMath.isPointInBox(new Vector3(1, 10, 10), new Vector3(0, 0, 0), new Vector3(1, 1, 1), new Quaternion(1, 1, 0, 0))).toBe(false)
    expect(AMath.isPointInBox(new Vector3(0.5, -0.5, 0.5), new Vector3(0, 0, 0), new Vector3(1, 1, 1), new Quaternion(Math.cos(Math.PI / 4), 1, 0, 0))).toBe(true)
})